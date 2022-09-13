# pylint: disable-msg=C0103, C0101
# -*- coding: ISO-8859-1 -*-

# Copyright (c) 2002 LOGILAB S.A. (Paris, FRANCE).
# http://www.logilab.fr/ -- mailto:contact@logilab.fr
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

"""Hidden Markov Models in Python
Implementation based on _A Tutorial on Hidden Markov Models and Selected
Applications in Speech Recognition_, by Lawrence Rabiner, IEEE, 1989,
 _Improved Estimation of Hidden Markov Model Parameters from Multiple 
Observation Sequences_, by Richard I. A. Davis, Brian C. Lovell, Terry Caelli,
2002, _Improved Ensemble Training for Hidden Markov Models Using Random 
Relative Node Permutations_, by Richard I. A. Davis and Brian C. Lovell, 2003.
This module is an heritage of the module HMM and uses numeric python multyarrays 
to improve performance and reduce memory usage 

Updated by MTD, 8/18/09:
Allows for multiple 'quality scores' (could be other things) for the observation
probabilities.  So instead of M * N for the observation_proba matrix, now it is
Q * M * N.  
In my case, with each observation set, I pass in the associate set of quality values.
These are an integer between 0 and 40.  
"""

from numpy import array, ones, zeros, add, cumsum, searchsorted, \
     product, dot, multiply, alltrue, log, equal, newaxis, \
     take, empty_like, allclose, where, sum
from numpy.random import random
import pickle

# Display log likelihood every DISPITER iterations while learning
DISPITER = 100

matrixproduct = dot

EPSILON = 1e-9
SMALLESTFLOAT = 1e-320
# These are the tolerances used with 'allclose()'
# For now we don't take the size of the model into
# account, but we really should. That's why you
# may want to tune this according to the size
# of your model.
# Since we use probability matrices the absolute
# tolerance will almost always come into play
# allclose(x,y) uses |x-y|<ATOL+RTOL*|y|
# because RTOL*|y| will be of the same order of magnitude
# as ATOL
# One more note, for some applications, the B matrix will be
# a lot larger than A. This is why we check A, and pi first
# and then B. This is also why you may want to tune tolerances
# so that alpha is less likely to trigger the test
alpha_ATOL = 1e-9
alpha_RTOL = 1e-6
beta_ATOL = 1e-8
beta_RTOL = 1e-5

## HMM Helper functions
## These functions are defined outside the HMM class
## because they should have the same prototype as their
## C counterpart

def _alpha_scaled(A, Bo, pi):
    """Internal method.
    Computes forward probabilities values, using a rescaling methods
    alpha_scaled[t,i]=Normalisation(P(O(1)...O(t),Q(t)=Si|model))
    Bo is the "slice" of the observation probability matrix corresponding
    to the observations (ie Bo=take(B,observation_indices)).
    For each t, c(t)=1./sum(alpha(t,i)), and C(t)=product(k=0..t,c(t))
    and alpha_scaled(t,i)=alpha(t,i)*C(t)
    The function returns: (alpha_scaled,C(t))
    """
    T = Bo.shape[0]
    N = A.shape[0]
    alpha_t = Bo[0] * pi                # (19)
    scaling_factors = zeros( T, float )
    scaling_factors[0] = 1./add.reduce(alpha_t)    
    alpha_scaled = zeros( (T, N), float)
    alpha_scaled[0] = alpha_t*scaling_factors[0]
    for i in range(1, T):
        alpha_t = dot(alpha_scaled[i-1], A)*Bo[i]  # (92a)        
        scaling_t = 1./add.reduce(alpha_t)
        scaling_factors[i] = scaling_t
        alpha_scaled[i] = alpha_t*scaling_t      # (92b)    
    return alpha_scaled, scaling_factors


def _beta_scaled( A, Bo, scale_factors ):
    """Computes backward probabilities
    beta(t,i)=P(O(t+1),...,O(T),Q(t)=Si|model)
    Or if scale_factors is not None:
    beta_scaled(t,i)=beta(t,i)*C(t) (From the result of _alpha_scaled)
    Bo is the same as in function _alpha
    """
    T, N = Bo.shape
    assert N == A.shape[0]
    scale_factors = scale_factors
    beta = zeros( (T, N), float )
    tmp = zeros( N, float )
    beta[-1] = ones( N, float ) * scale_factors[-1]         # (24)
    for t  in range( T-2, -1, -1 ):
        multiply( scale_factors[t], Bo[t+1], tmp )
        multiply( tmp, beta[t+1], tmp )
        beta[t] = matrixproduct( A, tmp )    # (25)
    return beta


def _ksi( A, Bo, alpha, beta ):
    """Compute ksi(t,i,j)=P(q_t=Si,q_(t+1)=Sj|model)"""
    N = A.shape[0]
    T = len(Bo)
    ksy = zeros( (T-1, N, N), float )
    tmp = Bo * beta
    for t in range(T-1):
        # This does transpose[alpha].(B[obs]*beta[t+1])
        # (where . is matrixproduct)
        ksit = ksy[t, :, :]
        multiply( A, tmp[t+1], ksit )
        multiply( ksit, alpha[t, :, newaxis], ksit )
        ksi_sum = add.reduce( ksit.flat )
        ksit /= ksi_sum
    return ksy

def _update_iter_B( gamma, qualList, obsIndices, B_bar ):
    """Updates the estimation of the observations probabilities.
    This function is used during the learning process."""
    # Contrary to the equations in the paper from rabiner
    # For B we sum over all the observations.
    # We cannot do this for A because it doesn't make sense
    # But for B we need to take into account the last symbol
    # in the chain (If we don't it leads to such things as
    # the probability of a fullstop at the end of a sentence
    # is zero!!)
    for i in range(len(obsIndices)):     # (110) numerateur
        B_bar[qualList[i], obsIndices[i]] += gamma[i]
        #B_bar[obsIndices[i]] += gamma[i]

def _correct_M( M, k, p ):
    """This function is a hack. It looks for states with 0 probabilities, and
    changes this probability to a uniform probability. This avoids divide by zero
    errors, and doesn't change the result of the algorithm.
    You can only have 0 probabilities if your observation matrix contains symbols
    that don't appear in your observations AND the initial state transition and
    observation probabilities are such that a state is reachable only if you observe
    those symbols.
    Parameters are:
    M the matrix
    k the axis along which we need a pdf
    p the value to replace with (usually 1/M.shape[k])
    """
    D = equal( add.reduce( M, k ), 0.0)
    if k == 1:
        for i in range(M.shape[0]):
            if D[i]:
                M[i, :] = p
    elif k == 0:
        for i in range(M.shape[1]):
            if D[i]:
                M[:, i] = p
    else:
        raise "Not Implemented"
    return M


def _normalize_B( B_bar, sigma_gamma_B ):
    """Internal function.
    Normalize the estimations of matrix A.
    Make sure we get rid of lines that contains only zeroes."""
    #print "Before _normalize_B"
    #print B_bar
    #print "-" * 60

    #sigma_gamma_B = 1./where( sigma_gamma_B, sigma_gamma_B, 1)
    #B_bar *= sigma_gamma_B    # (110)
    for i in range(len(B_bar)):
        total1 = B_bar[i][0][0] + B_bar[i][1][0]
        if total1 == 0:
            total1 = 1.0
        B_bar[i][0][0] /= total1
        B_bar[i][1][0] /= total1
        
        total2 = B_bar[i][0][1] + B_bar[i][1][1]
        if total2 == 0:
            total2 = 1.0
        B_bar[i][0][1] /= total2
        B_bar[i][1][1] /= total2
        
    #print "AFTER _normalize_B"
    #print B_bar
    #print "*"*60
    #print


## ----------------------------------------------------------------------

class HMM:
    """A Hidden Markov Model implementation
    Methods are provided for computing the probabitility of a sequence of
    observation, computing the most probable state transitions leading to
    a sequence of observations, as well as supervised and unsupervised
    training of the model and generation of sequences of observations
    (simulation).
    The notations (member variables and some method names), especially greek
    letter names are directly inspired by [Rabiner89] mentionned above.
    Comments in the source code mentionning a number are references to
    equations in the algorithm descriptions of that paper."""

    alpha_scaled = staticmethod(_alpha_scaled)
    beta_scaled = staticmethod(_beta_scaled)
    ksi = staticmethod(_ksi)
    update_iter_B = staticmethod(_update_iter_B)
    correct_M = staticmethod(_correct_M)
    normalize_B = staticmethod(_normalize_B)

    ORDER = "C"
    
    def __init__(self, state_list, observation_list,
                 num_quality_vals,
                 transition_proba = None,
                 observation_proba = None,
                 initial_state_proba = None):
        """Builds a new Hidden Markov Model
        state_list is the list of state symbols [q_0...q_(N-1)]
        observation_list is the list of observation symbols [v_0...v_(M-1)]
        transition_proba is the transition probability matrix
            [a_ij] a_ij = Pr(X_(t+1)=q_j|X_t=q_i)
        observation_proba is the observation probablility matrix
            [b_ik] b_ik = Pr(O_t=v_k|X_t=q_i)
        initial_state_proba is the initial state distribution
            [pi_i] pi_i = Pr(X_0=q_i)"""
        self.num_quality_vals = num_quality_vals
        self.N = len(state_list)
        self.M = len(observation_list)
        self.omega_X = state_list
        self.omega_O = observation_list
        if transition_proba is None:
            transition_proba = ones( (self.N, self.N), float) / self.N
        if observation_proba is None:
            observation_proba = ones( (num_quality_vals, self.M, self.N), float) / self.M
        if initial_state_proba is None:
            initial_state_proba = ones( (self.N,), float ) / self.N

        self.A = array(transition_proba, float, order=self.ORDER)
        self.B = array(observation_proba, float, order=self.ORDER)
        self.pi = array(initial_state_proba, float, order=self.ORDER)
        # dimensional assertions
        self.checkHMM()
        self.make_indexes()
        
    def make_indexes(self):
        """Creates the reverse table that maps states/observations names
        to their index in the probabilities array"""
        self.X_index = {}
        for i in range(self.N):
            self.X_index[self.omega_X[i]] = i
        self.O_index = {}
        for i in range(self.M):
            self.O_index[self.omega_O[i]] = i
            
    def saveHMM( self, f, saveState = None ):
        """Save the data for this class using cPickle.
        NOTE: don't use cPickle directly if your data uses
        too much memory. The pickle implementation of arrays
        just (well, not exactly) does a big binary copy of itself
        into a string and let pickle save the string object.
        So USE this function if your data becomes too big.
        As a side note, pickle will fail anyway because we have
        function objects as members of the HMM object. To use pickle
        you need to define __getattr__ __setattr__.
        """
        version = "HMM1.0"
        pickle.dump( version, f, 1 )
        pickle.dump( saveState, f, 1 )
        if saveState:
            pickle.dump( self.omega_X, f, 1 )
            pickle.dump( self.omega_O, f, 1 )
        pickle.dump( self.N, f, 1 )
        pickle.dump( self.M, f, 1 )
        pickle.dump( self.num_quality_vals, f, 1)
        pickle.dump( self.A, f, 1 )
        pickle.dump( self.pi, f, 1 )
        for i in range(self.num_quality_vals):
            for j in range(self.M):
                pickle.dump( self.B[i, j, :], f, 1 )
        
    def loadHMM( self, f ):
        """Use this function if you saved your data using
        saveHMM."""
        version = pickle.load(f, encoding='latin1')
        if version == "HMM1.0":
            saveState = pickle.load(f, encoding='latin1')
            if saveState:
                self.omega_X = pickle.load(f, encoding='latin1')
                self.omega_O = pickle.load(f, encoding='latin1')
            self.N = pickle.load(f, encoding='latin1')
            self.M = pickle.load(f, encoding='latin1')
            self.num_quality_vals = pickle.load(f, encoding='latin1')
            if saveState:
                self.make_indexes()
            self.A = pickle.load(f, encoding='latin1')
            self.pi = pickle.load(f, encoding='latin1')
            self.B = zeros( (self.num_quality_vals, self.M, self.N), float, self.ORDER )
            for i in range(self.num_quality_vals):
                for j in range(self.M):
                    x = pickle.load(f, encoding='latin1')
                    self.B[i, j, :] = x
        else:
            raise RuntimeError("File format not recognized")

    def checkHMM(self):
        """This function will asserts if the internal state of the class
        is inconsistent. (Checks that the matrices' sizes are correct and
        that they represents probabilities)."""
        assert self.A.shape == (self.N, self.N), \
               """transition_proba must be a N*N matrix, where N is 
                len(state_list)"""
        assert self.pi.shape == (self.N, ), \
               """transition_proba must be a N element vector,
               where N is len(state_list)"""
        reduced = add.reduce(self.A, 1) - 1
        assert (alltrue(reduced < EPSILON) and \
                alltrue(reduced > -EPSILON) and \
                alltrue(alltrue(self.A<=1.0)) and \
                alltrue(alltrue(self.A>=0.0))),\
                """transition_proba must be a probability matrix"""
        for oneB in self.B:
            assert oneB.shape == (self.M, self.N), \
               """observation_proba must be a Q*M*N matrix, where N is 
                len(state_list) and M is len(observation_list) and Q is
                the number of quality values.  This failed for %s""" % oneB
            reduced = add.reduce(oneB, 0) - 1
            assert (alltrue(reduced < EPSILON) and \
                alltrue(reduced > -EPSILON) and \
                alltrue(alltrue(self.B<=1.0)) and \
                alltrue(alltrue(self.B>=0.0))),\
                """each column of each quality's observation_proba must be a probability
                vector"""
        if len(self.pi)==0: # a zero length vector is reduced to a scalar
            return          # and makes the following test fail
        reduced = add.reduce(self.pi) - 1
        assert (reduced < EPSILON and reduced > -EPSILON and \
                alltrue(self.pi<=1.0) and \
                alltrue(self.pi>=0.0)), \
                """initial_state_proba must be a probability vector"""


    def dump(self):
        """Helper method for debugging"""
        print(self.getDumpStr())
        
        
    def getDumpStr(self):
        s = ""
        s += "=" * 80
        s += "\n%s\n" % (str(self))
        s += "States: %s\n" % self.N
        s += "Observations: %s\n" % self.M
        s += "-" * 80
        s += "\nState transition probabilities:\n"
        s += "\n%s\n" % (str(self.A))
        s += "-" * 80
        s += "\nObservation probabilities:\n"
        s += "%s\n" % (str(self.B))
        s += "-" * 80
        s += "\n"
        
        return s
        
    def dumpPretty(self):
        """Dumps in a single tab delimited line for easy import into excel for analysis"""
        vals = []
        for x in self.A:
            for y in x:
                vals.append(str(y))
        for oneB in self.B:
            for x in oneB:
                for y in x:
                    vals.append(str(y))
            
        print("|".join(vals))
        
    def __getinitargs__(self):
        """helper method for pickling"""
        return self.omega_X, self.omega_O, self.A, self.B, self.pi

    def set_random_proba(self):
        """Assigns random probability to all three matrices A, B and pi"""
        self.set_random_transition_proba()
        self.set_random_observation_proba()
        self.set_random_initial_proba()
        self.checkHMM()

    def reset_transition_proba(self):
        """This resets the state transition matrix to zero. Use it
        only if you want to use set_transition_proba on some coefficients."""
        multiply( self.A, 0.0, self.A )

    def reset_observation_proba(self):
        """This resets the observation matrix to zero. Use it
        only if you want to use set_observation_proba on some coefficients."""
        multiply( self.B, 0.0, self.B )

    def reset_initial_proba(self):
        """This resets the initial distribution matrix to zero. Use it
        only if you want to use set_initial_proba on some coefficients."""
        multiply( self.pi, 0.0, self.pi )

    def set_random_transition_proba(self):
        """set transition probability matrix to some random values"""
        self.A = random( self.A.shape )
        self.A /= add.reduce( self.A, 1 )[:, newaxis] # normalization

    def set_random_observation_proba(self):
        """set observation probability matrix to some random values"""
        self.B = random( self.B.shape )
        self.B /= add.reduce( self.B ) # normalization

    def set_random_initial_proba(self):
        """set initial state probability matrix to some random values"""
        self.pi = random( self.pi.shape )
        self.pi /= add.reduce( self.pi ) # normalization

    def set_transition_proba( self, state1, state2, value ):
        """set the probability of a transition form 'state1' to 'state2'"""
        self.A[ self.X_index[state1], self.X_index[state2] ] = value

    def set_observation_proba( self, quality, state, obs, value ):
        """set the probability of generating observation 'obs'
        when in state 'state'"""
        self.B[ quality ][self.O_index[obs], self.X_index[state] ] = value

    def set_initial_proba( self, state, value ):
        """set the probability of being initially in state 'state'"""
        self.pi[self.X_index[state]] = value

    def _get_observationIndices( self, observations ):
        """return observation indices"""
##        return [self.O_index[o] for o in observations]
        indices = zeros( len(observations), int )
        k = 0
        for o in observations:
            indices[k] = self.O_index[o]
            k += 1
        return indices
    
#    def simulate( self, length, show_hidden = False ):
#        """generates a random sequence of observations of given length
#        if show_hidden is true, returns a list of (state,observation)"""
#        cumA = cumsum( self.A, 1 )
#        cumB = cumsum( self.B, 0 )
#        r0 = random()
#        state = searchsorted( cumsum(self.pi, 0), r0)
#        seq = []
#        states = []
#        
#        for i in xrange(length):
#            states.append(state)
#            r1 = random()
#            symbol = self.omega_O[ searchsorted( cumB[:, state], r1 ) ]
#            if show_hidden:
#                seq.append( (self.omega_X[state], symbol) )
#            else:
#                seq.append(symbol)
#            r2 = random()
#            state = searchsorted( cumA[state, :], r2 )
#        return seq

    def analyze( self, observations , qualities):
        """use Viterbi algorithm to
        find the states corresponding to the observations.
        qualities is a list the length of the observations
        that contains the quality scores"""
        B = self.B
        A = self.A
        T = len(observations)
        N = self.N
        Omega_X = self.omega_X
        obs = self._get_observationIndices(observations)
        q = [int(x) for x in qualities]
        # initialisation
        delta = zeros( N, float )
        tmp = zeros( N, float )
        delta = B[q[0], obs[0]] * self.pi    # (32a)
        delta_t = zeros( N, float )
        psi = zeros( (T, N), int )       # (32b)
        # recursion
        for t in range(1, T):
            O_t = obs[t]
            q_t = q[t]
            for j in range(N):
                multiply( delta, A[:, j], tmp )
                idx = psi[t, j] = tmp.argmax()        # (33b)
                delta_t[j] = tmp[idx] * B[q_t, O_t, j]  # (33a)
            delta, delta_t = delta_t, delta

        # reconstruction
        i_star = [delta.argmax()]                         # (34b)
        for psi_t in psi[-1:0:-1]:
            i_star.append( psi_t[i_star[-1]] )                 # (35)
        trajectory = [Omega_X[i] for i in i_star]
        trajectory.reverse() # put time back in the right direction
        return trajectory
        
#    def analyze_log( self, observations ):
#        """use a modified Viterbi algorithm (using log P) to
#        find the states corresponding to the observations."""
###      Since we use log(), we replace log(0) by the variable
###      hmm.SMALLESTFLOAT, change it if you find it not small enough. 
#        B = self.B
#        A = self.A
#        T = len(observations)
#        N = self.N
#        M = self.M
#        Omega_X = self.omega_X
#        obs = self._get_observationIndices(observations)
#        k = equal( A, 0.0 ) * SMALLESTFLOAT
#        logA = log( A + k )
#        logB = zeros( (M, N), float)
#        for i in xrange(M):
#            t = B[i, :]
#            k = equal( t, 0.0 ) * SMALLESTFLOAT
#            logB[i] = log( k + t )
#        # initialisation
#        psi = zeros( N, float )
#        psi_t = zeros( N, float )
#        logPi = log( self.pi + equal( self.pi, 0.0 ) * SMALLESTFLOAT )
#        add( logB[obs[0]], logPi, psi) # (105a)
#        Q = zeros( (T, N), int )
#        # recursion
#        tmp = zeros( N, float )
#        for t in xrange( 1, T ):
#            O_t = obs[t]
#            for j in xrange(N):
#                tmp = psi + logA[:, j]
#                idx = Q[t, j] = tmp.argmax()
#                psi_t[j] = tmp[idx] + logB[O_t, j] # (105b)
#            psi, psi_t = psi_t, psi
#
#        # reconstruction
#        q_star = [psi.argmax()]                         # (105c)
#        for q_t in Q[-1:0:-1]:
#            q_star.append( q_t[q_star[-1]] )                 # (35)
#        trajectory = [Omega_X[i] for i in q_star]
#        trajectory.reverse() # put time back in the right direction
#        return trajectory

    def log_likelihood( self, observations, qualities, trajectory ):
        """return log_likelihood"""
        obs = self._get_observationIndices(observations)
        states = [ self.X_index[s] for s in trajectory ]
        q = [int(x) for x in qualities]
        res = 0
        N = self.N
        M = self.M
        logB = zeros( (M, N), float)
        for i in range(M):
            t = self.B[q[i], i, :]
            k = equal(t, 0.0) * SMALLESTFLOAT
            logB[i] = log(k + t)
        for o, s in zip( obs, states ):
            res += logB[o, s]
        return res

    def _mask(self):
        """forgive round errors"""
        mask = (1 < self.A) & (self.A  <= 1+EPSILON)
        self.A[mask] = 1.0
        mask = (1 < self.B) & (self.B  <= 1+EPSILON)
        self.B[mask] = 1.0
        mask = (1 < self.pi) & (self.pi  <= 1+EPSILON)
        self.pi[mask] = 1.0

#    def learn( self, observations, maxiter = 1000, impr=1 ):
#        """Train the model according to one sequence of observations"""
#        obs = self._get_observationIndices(observations)
#        iter = self._baum_welch( obs, maxiter, impr )
#        return iter
#
#    def _get_observations( self, obsIndices ):
#        """Extract the lines of the observations probability matrix 
#        corresponding to the actual observations."""
#        return take( self.B, obsIndices, 0 )
#
    def _likelihood( self, scale_factors ):
        """This function computes the log likelihood
        of the training set using the precomputed
        alpha probabilities (sum(k=0..N,alpha(T,k)).
        It should increase during the learning process."""
        t = where( scale_factors==0.0, SMALLESTFLOAT, scale_factors )
        return -add.reduce( log(t) )

    def multiple_learn(self, m_observations, m_quals,
                       maxiter = 1000, impr=1 ):
        """Uses Baum-Welch algorithm to learn the probabilities on multiple
        observations sequences
        """
        # remove empty lists
        m_observations = [x for x in m_observations if x]
        setO =  set()   # set of observations        
        K = len( m_observations )
        learning_curve = []
        sigma_gamma_A = zeros( (self.N, ), float, order=self.ORDER )
        sigma_gamma_B = zeros( (self.N, ), float, order=self.ORDER )
        A_bar  = zeros( (self.N, self.N), float, order=self.ORDER )
        B_bar  = zeros( (self.num_quality_vals, self.M, self.N), float, order=self.ORDER )
        pi_bar = zeros( self.N, float, order=self.ORDER )
        if DISPITER == 0:
            dispiter = maxiter
        else:
            dispiter = DISPITER
        obs_list = []
        for k in range(K):
            observations = m_observations[k]
            obsIndices = self._get_observationIndices(observations)
            obs_list.append( obsIndices )
            setO = setO | set(obsIndices)  # add new elements observed
        for iter in range( 1, maxiter + 1 ):
            # mtd
            #if impr and iter % 10 == 0:
            #    print iter
            #    self.dump()
            total_likelihood = 0
            for k in range(K):
                obsIndices = [int(j) for j in obs_list[k]]
                qualList = [int(j) for j in m_quals[k]]
                #Bo = take(self.B, obsIndices, 0)
                # Bo is an array with each state's proba depending on the observations.  
                # For now, build by hand from both quality and state info
                Bo = zeros( (len(obsIndices), self.N), float)
                for i in range(len(obsIndices)):
                    Bo[i] = self.B[qualList[i], obsIndices[i]]
                alpha, scale_factors = self.alpha_scaled( self.A, Bo, self.pi )
                beta  = self.beta_scaled( self.A, Bo, scale_factors )
                ksy   = self.ksi( self.A, Bo, alpha, beta )
                gamma = self._gamma( alpha, beta, scale_factors )
                pi_bar += gamma[0]
                self._update_iter_gamma( gamma, sigma_gamma_A, sigma_gamma_B )
                self._update_iter_A( ksy, A_bar )
                self.update_iter_B( gamma, qualList, obsIndices, B_bar )
                total_likelihood += self._likelihood( scale_factors )
                
            #end for k in range(K)
            self._normalize_iter_A( A_bar, sigma_gamma_A )
            self.normalize_B( B_bar, sigma_gamma_B )
            pi_bar /= K
            # Correct A_bar and B_bar in case 0 probabilities slipped in
            self.correct_M( A_bar, 1, 1. / self.N )
            # TODO: add this back in
            #self.correct_M( B_bar, 0, 1. / self.M )
            learning_curve.append( total_likelihood )
            if (iter % dispiter) == 0:
                if impr:
                    print("Iter ", iter, " log=", total_likelihood)
            if self._stop_condition( A_bar, pi_bar, B_bar ):
                if impr:
                    print('Converged in %d iterations' % iter)
                break
            self.A, A_bar   = A_bar, self.A
            self.B, B_bar   = B_bar, self.B
            self.pi, pi_bar = pi_bar, self.pi
            A_bar.fill(0)
            B_bar.fill(0)
            pi_bar.fill(0)
            sigma_gamma_A.fill(0)
            sigma_gamma_B.fill(0)
        else:
            if impr:
                print("The Baum-Welch algorithm did not converge in", end=' ')
                print(" %d iterations" % maxiter)
        self._mask()
        # Correct B in case 0 probabilities slipped in
        setO = set(range(self.M)) - setO
        while setO:
            e = setO.pop()
            self.B[e] = 0
        return iter, learning_curve

#    def _baum_welch( self, obsIndices, maxiter=1000, impr=1 ):
#        """Uses Baum-Welch algorithm to learn the probabilities
#        Scaling on the forward and backward values is automatically
#        performed when numerical problems (underflow) are encountered.
#        Each iteration prints a dot on stderr, or a star if scaling was
#        applied"""
#        B  = self.B
#        A  = self.A
#        pi = self.pi
#        learning_curve = []
#        if DISPITER == 0:
#            dispiter = maxiter
#        else:
#            dispiter = DISPITER
#        Bo = take( B, obsIndices, 0 )
#        for iter in xrange( 1, maxiter + 1 ):
#            alpha, scale_factors = self.alpha_scaled( A, Bo, pi )
#            beta = self.beta_scaled( self.A, Bo, scale_factors )
#            ksy = self.ksi( self.A, Bo, alpha, beta )
#            gamma = self._gamma( alpha, beta, scale_factors )
#            A_bar, B_bar, pi_bar = self._final_step( gamma, ksy, obsIndices )
#            learning_curve.append( self._likelihood(scale_factors) )
##            if impr:
##                if (iter % dispiter) == 0:
##                    print "Iter ", iter, " log=", learning_curve[-1]
#            if self._stop_condition( A_bar, pi_bar, B_bar):
##                if impr:
##                    print 'Converged in %d iterations' % iter
#                break
#            else:
#                self.A = A = A_bar
#                self.B = B = B_bar
#                self.pi = pi = pi_bar
#                self._mask()
#        else:
#            if impr:
#                print "The Baum-Welch algorithm did not converge",
#                print " in %d iterations" % maxiter
#        return iter, learning_curve
#
    def _update_iter_gamma( self, gamma, sigma_gamma_A, sigma_gamma_B ):
        """update iter gamma"""
        sigma_gamma_kA = add.reduce(gamma[:-1])
        sigma_gamma_A += sigma_gamma_kA       # (109) et (110) denominateur
        sigma_gamma_B += sigma_gamma_kA + gamma[-1]

    def _update_iter_A( self, ksy, A_bar ):
        """update iter A"""
        A_bar_k = add.reduce( ksy )
        add( A_bar, A_bar_k, A_bar )           # (109) numerateur
    
    def _normalize_iter_A( self, A_bar, sigma_gamma_A ):
        """Internal function.
        Normalize the estimations of matrix A.
        Make sure we get rid of lines that contains only zeroes."""
        # replace 0 with 1 to avoid div error
        # it doesn't matter if it is one or anything else since
        # sigma_gamma(i)=0 implies A(i,:)=0 and B(i,:)=0
        sigma_gamma_A = 1. / where( sigma_gamma_A, sigma_gamma_A, 1 )
        A_bar *= sigma_gamma_A[:, newaxis]    # (109)
#
#    def _final_step( self, gamma, ksy, obsIndices ):
#        """Compute the new model, using gamma and ksi"""
#        sigma_gamma_A = add.reduce(gamma[:-1])
#        sigma_gamma_B = add.reduce(gamma)
#        for i in range(len(sigma_gamma_B)):
#            if sigma_gamma_B[i] < EPSILON:
#                sigma_gamma_B[i] = 1
#        for i in range(len(sigma_gamma_A)):
#            if sigma_gamma_A[i] < EPSILON:
#                sigma_gamma_A[i] = 1
#        ## Compute new PI
#        pi_bar = gamma[0]                       # (40a)
#        ## Compute new A
#        A_bar  = add.reduce(ksy)
#        A_bar /= sigma_gamma_A[:, newaxis] # (40b)       
#        ## Compute new B
#        B_bar = zeros( (self.M, self.N), float )
#        for i in xrange( len(obsIndices) ):
#            B_bar[obsIndices[i]] += gamma[i] 
#        B_bar /= sigma_gamma_B
#        return A_bar, B_bar, pi_bar
       
    def _stop_condition( self, A_bar, pi_bar, B_bar ):
        """Returns true if the difference between the estimated model
        and the current model is small enough that we can stop the
        learning process"""
        return (allclose( self.A, A_bar, alpha_RTOL, alpha_ATOL) and 
               allclose( self.pi, pi_bar, alpha_RTOL, alpha_ATOL) and 
               allclose( self.B, B_bar, beta_RTOL, beta_ATOL))
    
    def _gamma(self, alpha, beta, scaling_factors ):
        """Compute gamma(t,i)=P(q_t=Si|model)"""
        g = alpha * beta / scaling_factors[:, newaxis]
        return g

#    def normalize(self, P = None):
#        """This can be used after a learning pass to
#        reorder the states so that the s_i -> s_i transition
#        probability are ordered s(0,0)>s(1,1) ... > s(n,n)
#
#        the permutation of states can be passed as a parameter
#        """
#        if P is None:
#            P = self.A.diagonal().argsort()
#            P = P[::-1]
#        A = empty_like(self.A)
#        PI = empty_like(self.pi)
#        B = empty_like(self.B)
#        N = A.shape[0]
#        for i in xrange(N):
#            pi = P[i]
#            for j in xrange(N):
#                pj = P[j]
#                A[i, j] = self.A[pi, pj]
#            B[:, i] = self.B[:, pi]
#            PI[i] = self.pi[pi]
#        return A, B, PI
#
#    def _weighting_factor_Pall(self, setObs):
#        """compute Wk = P(setObservations | lambda_k) """
#        P = 1
#        for obs in setObs:
#            Tk = len(obs)
#            obsIndices = self._get_observationIndices(obs)
#            Bo = take(self.B, obsIndices, 0)
#            null = 0
#            for i in range(Tk):
#                null = null or (allclose(Bo[i], zeros([self.N])))
#            if null:
#                P = 0
#            else:
#                alpha_s, scalingFactor = self.alpha_scaled(self.A, Bo, self.pi)
#                alpha = alpha_s[Tk-1] / product(scalingFactor, 0) 
#                P *= add.reduce(alpha)
#        return P
#
#    def _weighting_factor_Pk(self, observation):
#        """compute Wk = P(Observation_k | lambda_k) """
#        Tk = len(observation)
#        obsIndices = self._get_observationIndices(observation)
#        Bo = take(self.B, obsIndices, 0)
#        alpha_s, scalingFactor = self.alpha_scaled(self.A, Bo, self.pi)
#        alpha = alpha_s[Tk-1] / product(scalingFactor, 0)
#        return add.reduce(alpha)
#
#    def ensemble_averaging(self, setObservations, weighting_factor="unit", 
#                            maxiter=1000, impr=1):
#        """Uses ensemble averaging method to learn the probabilities on multiple
#        observations sequences"""
#        N = self.N
#        W = 0
#        self._mask()
#        hmmk = HMM(self.omega_X, self.omega_O, self.A, self.B, self.pi)
#        A_bar = zeros( (N, N))
#        B_bar = zeros( (self.M, N))
#        pi_bar = zeros(N)
#        for obs in setObservations:
#            hmmk.A = self.A
#            hmmk.B = self.B
#            hmmk.pi = self.pi
#            obsIndices = self._get_observationIndices(obs)
#            hmmk._baum_welch(obsIndices, maxiter, impr)
#            if weighting_factor == "Pall":
#                Wk = hmmk._weighting_factor_Pall(setObservations)
#            elif weighting_factor == "Pk":
#                Wk = hmmk._weighting_factor_Pk(obs)
#            else:
#                Wk = 1
#            A_bar = A_bar + Wk * hmmk.A
#            B_bar = B_bar + Wk * hmmk.B
#            pi_bar = pi_bar + Wk * hmmk.pi
#            W = W + Wk
#        if W == 0:
#            W = 1
#            print "The ensemble averaging method did not converge"
#        else:
#            self.A = A_bar / W
#            self.B = B_bar / W
#            self.pi = pi_bar / W
#            self._mask()
#            
#def testLearn():
#    test = HMM(['a','b'],['s1','s2','s3'])
#    test.set_random_proba()
#    
#    #test = HMM(['a','b'],['s1','s2','s3'],
#    #           array([[.3,.5],[.7,.5]]),
#    #           array([[.45,.1],[.45,.45],[.1,.45]]),
#    #           array([.9,.1]))
#    
#    #observe = array( [ ['s1', 's1', 's2', 's3'], ['s1', 's2', 's2', 's3'], ['s1', 's1', 's2', 's3'], 
#    #                   ['s1', 's1', 's2', 's3'], ['s1', 's2', 's2', 's3'] ] )
#    
#    observe = array( [ ['s1', 's1', 's1', 's1'], ['s2', 's2', 's1', 's3'], ['s3', 's1', 's2', 's3'], 
#                       ['s3', 's1', 's2', 's1'], ['s3', 's2', 's1', 's3'] ] )
#    
#    print "BEFORE:"
#    print "A:"
#    print test.A
#    print "B"
#    print test.B
#    
#    test.multiple_learn(observe, 100)
#    
#
#    print "AFTER:"
#    print "A:"
#    print test.A
#    print "B"
#    print test.B

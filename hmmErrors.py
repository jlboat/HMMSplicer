#    Copyright 2009,2010 Michelle Dimon
#
#    This file is part of HMMSplicer.
#
#    HMMSplicer is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HMMSplicer is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HMMSplicer.  If not, see <http://www.gnu.org/licenses/>.
#
class InvalidFileFormatException (Exception):
    pass

class CommandLineException (Exception):
    pass

class UnexpectedException(Exception):
    """Used for errors that 'should not happen'.  For example, if you've checked
    that the length is greater than some value and then you check later and the 
    length is shorter, you'd throw this exception."""
    pass

class InvalidInputException(Exception):
    """Raised when a function has been given invalid input.  For example, if a variable
    was expected and None was passed in."""
    pass

class InvalidQuality(Exception):
    """Raised when a function has been an invalid quality value."""
    pass

class InvalidFastq(Exception):
    """Raised when the input fastq format is invalid."""
    pass

class InvalidPythonVersion(Exception):
    """Raised when an invalid python version is used, for example Python 3.0.
    Also raised if the user tries to run multiple processors using Python 2.5.
    """
    pass
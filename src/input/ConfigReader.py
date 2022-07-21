
from eMolFrag2.src.utilities import logging, constants
from pathlib import Path
    
def grabCommands(contents):
  """
      Takes the string contents of a configuration file
      and grabs its command line arguments.
      If there are comments, ignore anything to the right of it
      
      @output: a list of flags and arguments
  """
  retString = ""
  position = 0

  #Read line by line, ignore comments, concatenate the remaining tokens into one string
  for line in contents:
      position = line.find("#")
      if (position >= 0):
          retString += line[:position]
      else:
          retString += line

  return retString.split(" ")

def readConfig(config_file, parser):
  """
      Reads a config file and parses arguments
      If the file is empty, throw an error
      
      @output: parsed arguments from argparser
  """
  config_file = Path(config_file)
  if config_file.suffix != constants.EMF_FORMAT_EXT:
    logging.logger.error(f"Configuration files must have the {constants.EMF_FORMAT_EXT} extension")
    return None

  with open (config_file) as f:
    lines = f.readlines()
  
  if (len(lines) <= 0):
    logging.logger.error(f"{config_file} is empty")
    return None

  #Grab the commands and then parse them
  cmdList = grabCommands(lines)

  args = parser.parse_args(cmdList)

  if args.i is None or args.o is None:
    f'Command-line arguments failed to parse; execution of eMolFrag will stop.'
    return None

  return args
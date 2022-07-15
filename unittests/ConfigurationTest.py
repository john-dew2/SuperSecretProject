import sys
from pathlib import Path
from eMolFrag2.unittests import utilities
from eMolFrag2.src.input import Options, Configuration

usr_dir = Path.cwd()
config_files = usr_dir.joinpath("eMolFrag2/unittests/data/config-files")

dataPath = usr_dir.joinpath("eMolFrag2/unittests/data/db-files")
mol2 = dataPath.joinpath("mol2")
smi = dataPath.joinpath("smi")
sdf = dataPath.joinpath("sdf")

def runReadCommandLine(arguments, expec_result):
    initializer = Options.Options()
    initializer = Configuration.readCommandLine(initializer, arguments)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadCommandLineTests():
    
    #Normal
    runReadCommandLine("-m eMolFrag2.src.eMolFrag -i {usr_dir.joinpath('eMolFrag2/unittests/data/db-files/mol2')} -o output/ -p 16 -c 1".split(" "), True)
    
    #Empty
    runReadCommandLine(None, False)
    
    
def runReadConfigurationFile(config_file, expec_result):
  assert Configuration.readConfigurationFile(config_file) == expec_result

def runReadConfigurationFileTests():
    #Files with a comment
    runReadConfigurationFile(config_files.joinpath("comment1.txt"), ['\n!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile(config_files.joinpath("comment2.txt"), ['!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1', ''])
    runReadConfigurationFile(config_files.joinpath("comment3.txt"), ['!python',  '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1\n\n'])
    runReadConfigurationFile(config_files.joinpath("comment4.txt"), [''])

    #Normal Files
    runReadConfigurationFile(config_files.joinpath("normal1.txt"),  ['!python',  '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '16', '-c', '1'])
    runReadConfigurationFile(config_files.joinpath("normal2.txt"),  ['!python', '-m', 'eMolFragTEMP.src.eMolFrag', '-i', '/content/eMolFragTEMP/test/mol2-test', '-o', 'output/', '-p', '0', '-c', '0'])

    #Empty File
    runReadConfigurationFile(config_files.joinpath("empty1.txt"), [])

def runReadConfigurationInput(arguments, expec_result):
    initializer = Options.Options()
    initializer = Configuration.readConfigurationInput(initializer, arguments)
    if (initializer != None):
      result = True
    else:
      result = False
    assert result == expec_result

def runReadConfigurationInputTests():
    #Normal processing
    runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -i {mol2} -o output/ -p 16 -c 1".split(" "), True)
    runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -i {smi} -o output/".split(" "), True)
    
    #No input or output
    runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -o output/ -p 16 -c 1".split(" "), False)
    runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -i {smi} -p 16 -c 1".split(" "), False)
    runReadConfigurationInput(f"-m eMolFrag2.src.eMolFrag -p 16 -c 1".split(" "), False)
    
    #Not enough arguments 
    runReadConfigurationInput(f"!python".split(" "), False)
    
    #Files only
        #Files only
    runReadConfigurationInput(f"-m {config_files.joinpath('comment1.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {config_files.joinpath('comment2.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {config_files.joinpath('comment3.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {config_files.joinpath('comment4.txt')}".split(" "), False)
    runReadConfigurationInput(f"-m {config_files.joinpath('empty1.txt')}".split(" "), False)
    runReadConfigurationInput(f"-m {config_files.joinpath('normal1.txt')}".split(" "), True)
    runReadConfigurationInput(f"-m {config_files.joinpath('normal2.txt')}".split(" "), True)
#
#
# Unit test initiation functionality
#
#
def runtest1(func):

    try:
       func()
       return True

    except:
       return False


def runtest(test_name, test_func, successful, failed):
    if runtest1(test_func):
        successful.append(test_name)
    else:
        failed.append(test_name)

def runtests():
    printlevel = 1
    utilities.emit(printlevel, f"Executing {__name__} unit tests.")
    
    #
    # Define all tests as a Dictionary: {str-name, <function-to-execute>}
    #
    tests = {"readCommandLine": runReadCommandLineTests, "readConfigurationFile": runReadConfigurationFileTests, "readConfigurationInput": runReadConfigurationInputTests}

    #
    # Run
    #
    successful = []
    failed = []
    for test_name, test_func in tests.items():
        runtest(test_name, test_func, successful, failed)

    # 
    # Report
    #
    if not failed:        
        utilities.emit(printlevel, f'{__name__} unit tests are successful.')

    else:
        for test in failed:
            utilities.emit(printlevel+1, f'Failed {test}.')

if __name__ == "__main__":
    runtests()
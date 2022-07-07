#Options.py


#input_path
#output_path
#parallel_cores_used
#outut_type:
    # if 0 then full_process
    # if 1 then chop_only
    # if 2 then chop_and_remove
#output_format:
    #if 0 then traditional_format
    #if 1 then different_files and remove_log_files
    #if 2 then different_files and !remove_log_files
    
class Options:
    def __init__(self):
        self.INPUT_PATH = ""
        self.OUTPUT_PATH = ""
        self.PARALLEL_CORES_USED = 0
        self.FULL_PROCESS = False
        self.CHOP_ONLY = False
        self.CHOP_AND_REMOVE = False
        #self.TRADITIONAL_FORMAT  = False
        self.DIFFERENT_FILES  = False
        self.REMOVE_LOG_FILES  = False

    #prints the current state of the options
    def printOptions(self):
        print("INPUT_PATH: ", self.INPUT_PATH, "\n")
        print("OUTPUT_PATH: ", self.OUTPUT_PATH, "\n")
        print("PARALLEL_CORES_USED: ", self.PARALLEL_CORES_USED, "\n")
        print("FULL_PROCESS: ", self.FULL_PROCESS, "\n")
        print("CHOP_ONLY: ", self.CHOP_ONLY, "\n")
        print("CHOP_AND_REMOVE: ", self.CHOP_AND_REMOVE, "\n")
        print("TRADITIONAL_FORMAT: ", self.TRADITIONAL_FORMAT, "\n")
        print("DIFFERENT_FILES: ", self.DIFFERENT_FILES, "\n")
        print("REMOVE_LOG_FILES: ", self.REMOVE_LOG_FILES, "\n")
    
    #prints when there's an error with inputs pertaining to -p -m or -c    
    def paramErr(self, argType, option, lower, upper):
        print("[Error] ", argType, " expects values between ", lower, "-", upper, " inclusive, your input: ", option)
       
    #takes in the arguments from the command line and updates each command's repsective option
    def handleOption(self, argType, option):
        if (argType == "-i"):   
            self.INPUT_PATH = option
            
        if (argType == "-o"):
            self.OUTPUT_PATH = option
            
        if (argType == "-p"):
            option = int(option)
            #throw error
            if ((option > 16) or (option < 0)):
                self.paramErr(argType, option, 0, 16)
            
            self.PARALLEL_CORES_USED = option
            
        if (argType == "-m"):
            option = int(option)
            #throw error
            if ((option > 2) or (option < 0)):
                self.paramErr(argType, option, 0, 2)
                
            if (option == 0):
                self.FULL_PROCESS = True
            if (option == 1):
                self.CHOP_ONLY = True
            if (option == 2):
                self.CHOP_AND_REMOVE = True
                
        if (argType == "-c"):
            option = int(option)
            #throw error
            if ((option > 2) or (option < 0)):
                self.paramErr(argType, option, 0, 2)
            
            if (option == 0):
                self.TRADITIONAL_FORMAT = True
            if (option == 1):
                self.DIFFERENT_FILES = True
                self.REMOVE_LOG_FILES = True
            if (option == 2):
                self.DIFFERENT_FILES = True
                self.REMOVE_LOG_FILES = False
                
    def parseCommandLine(self, arguments):
        argTypes = ["-i","-o","-p","-m","-c"] 

        #parse through each argument and send them into the handler
        for i in range(len(arguments)):
            if (arguments[i] in argTypes):
                self.handleOption(arguments[i], arguments[i + 1])
              
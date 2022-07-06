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
    
INPUT_PATH = ""
OUTPUT_PATH = ""
PARALLEL_CORES_USED = 0
FULL_PROCESS = False
CHOP_ONLY = False
CHOP_AND_REMOVE = False
TRADITIONAL_FORMAT  = False
DIFFERENT_FILES  = False
REMOVE_LOG_FILES  = False
    
def paramErr(argType, option, lower, upper)
    print("[Error] ", argType, " expects values between ", lower, "-", upper, " inclusive, your input: ", option)

def parseCommandLine(arguments):
    argTypes = ["-i","-o","-p","-m","-c"] 

    #parse through each argument and send them into the handler
    for i in range(len(arguments):
        if (argument in argTpyes):
            handleOption(arguments[i], arguments[i + 1])
    
def handleOption(argType, option)
    if (argType == "-i")
        INPUT_PATH = option
        
    if (argType == "-o")
        OUTPUT_PATH = option
        
    if (argType == "-p")
        option = (int)option
        
        if (option > 16) or (option < 0)):
            paramErr(argType, option, 0, 16)
        
        PARALLEL_CORES_USED = option
        
    if (argType == "-m")
        option = (int)option
        if (option > 2) or (option < 0)):
            paramErr(argType, option, 0, 2)
            
        if (option == 0)
            FULL_PROCESS = True
        if (option == 1)
            CHOP_ONLY = True
        if (option == 2)
            CHOP_AND_REMOVE = True
            
    if (argType == "-c")
        option = (int)option
        if (option > 2) or (option < 0)):
            paramErr(argType, option, 0, 2)
        
        if (option == 0)
            TRADITIONAL_FORMAT = True
        if (option == 1)
            DIFFERENT_FILES = True
            REMOVE_LOG_FILES = True
        if (option == 2)
            DIFFERENT_FILES = True
            REMOVE_LOG_FILES = False
    
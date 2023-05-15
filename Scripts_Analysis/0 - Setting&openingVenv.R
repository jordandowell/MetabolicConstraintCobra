###MAKING virtual environmetn
###
library("reticulate")
#virtualenv_create(envname  = "./MetabolismENV")
Sys.setenv(RETICULATE_PYTHON = "./MetabolismENV/bin/python3")

#activate envirnoment
use_virtualenv("./MetabolismENV")
py_config()






###MAKING virtual environmetn
###
library("reticulate")
#virtualenv_create(envname  = "./MetabolismENV")
Sys.setenv(RETICULATE_PYTHON = "./MetabolismENV/bin/python3")

#activate envirnoment
use_virtualenv("./MetabolismENV")
py_config()




<<<<<<< HEAD
=======

>>>>>>> 6140148fffe1d9cbb4038f7c24a3e0065bb5754c

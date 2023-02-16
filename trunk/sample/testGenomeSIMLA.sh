#!/bin/bash

###############################################################################
#                 genomeSIMLA example test commands                           #
###############################################################################

#module load genomeSIMLA/1.2.0

# generate dataset
genomeSIMLA metaDimSimTest.sim

# analyze based on generated dataset
genomeSIMLA metaDimSimTest.sim -l 450

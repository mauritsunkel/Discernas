
# Debugging and inspecting source toolkit
browser() # enter debugging where this line of code is added to the script
## can step into function, and functions within functions, also for generics!
recover() # calls browser() for a specific call-level
debug(func.name) # set debugging flag/breakpoint on first line on function --> to enter browser() at that line
debugonce(func.name) # debug() but only runs 1 time for repeated function calls
undebug(func.name) # remove debugging flag/breakpoint status from function
trace(func.name, at = n) # where n is line number of where you want to start debugging
trace(func.name, edit = T) # get pop-up window with editable function code, Q: DOES THIS CODE GET SAVED?
traceback() # shows traceback of last error/ran code, can also be seen after error

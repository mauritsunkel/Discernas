
beep <- function(n = 5) {
  for(i in seq(n)){
    system("rundll32 user32.dll, MessageBeep -1")
    Sys.sleep(.5)
  }
  print(paste0("beeped at: ", format(Sys.time(), "%F %H-%M-%S")))
}

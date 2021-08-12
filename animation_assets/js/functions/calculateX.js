const getFrequency = (file, mode) => {
  return float(file[mode]);
}
  
const getVibrationMode = (file, mode) => {
  return float(file[mode].split(' '))
}
  
  
const calculateX = (vibrationMode, frequency, time) => {
  x = []
  for(i = 0; i < vibrationMode.length; i++)
    x[i] = vibrationMode[i]*cos(frequency*time)
  return x 
}
  


 New command line options added to starge to v2

 NOTE in long options, fcommand line
 format is like "stargen --gasdust 50". not "stargen --gasdust50"
 not parameter name and paremeter together, not as like in sort options -
 "stargen -s1"

 option default_value
 
 stargen
  --gasdust 50.0
  --density 1.0
  --alpha 5.0
  --gamma 3.0
  --beta 1.2e-5
  --disk type 0 can be 1 , 2: 1
  --nearest 0.3
  --farthest 50.0
  --diskradius 200
  --cloud_eccentricity 0.2
  
  --migrate 1.0
  
  and untested, not OK:
  
  [--filterhill]
  [--filterperiod]    


--disk_type

0 normal , classic Dole-style disk


1 exponent distribution disk

 --disk_type 1 use with alpha
 like this
 
 stargen --disk-type 1 --alpha 4.0
 
 2: normal distribution disk, often "ring"
 
 alpha is here surface density exponent
 
  --disk_type 2 use with alpha ang gamma
  
   stargen --disk-type 2 --alpha 5.0 --gamma 1.0
   
   where alpha is standard deviation, "sigma"
   and gamma is mean location "mu"
   
   
 
 

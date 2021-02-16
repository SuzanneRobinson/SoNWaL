##Suehrcke solar radiation calculation
##K_clear values from https://www.sciencedirect.com/science/article/pii/S0960148114002857#tbl1


N=8 #day length in hours
n =20#day of the year, starting jan 1st as 1
#delta is the declination of the sun (180/pi to convert from degrees to radians)
delta=(23.45*sin(360*(284+n/365)))*180/pi
lat = 56.459*(180/pi)

K_clear = 0.579 # see paper for monthly values 

H_h =6 #monthly average of daily horizontal surface irradiation. (sunshine hours)

h_ss = acos(-tan(lat)*tan(delta))

I_sc = 1367 #solar constant and is equal to 1367 Wm−2


I_0 = I_sc*(1+0.033*cos(360/365)*n)


#Daily horizontal extra-terrestrial solar irradiation
H_o = (3600*24/pi)*I_0*(pi*h_ss/180)*(sin(lat)*sin(delta)+cos(lat)*cos(delta)*sin(h_ss))


#The monthly average daily clearness index
K = H_h/H_o

f_clear = (K/K_clear)^2

##Suehrcke solar radiation calculation
##K_clear values from https://www.sciencedirect.com/science/article/pii/S0960148114002857#tbl1


N=12 #day length in hours
n =20#day of the year, starting jan 1st as 1
S_h =6 #observed monthly average of daily horizontal surface irradiation. (sunshine hours)

#delta is the declination of the sun (pi/180 to convert from degrees to radians)
delta=(23.45*(sin(pi/180*(360*(284+n/365)))))
lat = (pi/180)*56.459

K_clear = 0.579 # see paper for monthly values 

#f_clear equiv to sunshine fraction S
S = S_h/N


h_ss = acos(pi/180*(-tan(lat)*tan(delta)))

I_sc = 1367 #solar constant and is equal to 1367 Wm−2


I_0 = I_sc*(1+0.033*cos(pi/180*((360/365)*n)))


#Daily horizontal extra-terrestrial solar irradiation
H_o = ((3600*24/pi)*I_0)*((pi*h_ss)/180)*((sin(lat)*sin(delta)+cos(lat)*cos(delta)*sin(h_ss*pi/180)))


#The monthly average daily clearness index
K = sqrt(S*K_clear)

H_h=K*H_o

H_h/(S_h*60*1000)






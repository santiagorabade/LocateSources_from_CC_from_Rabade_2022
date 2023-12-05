function [dist] = get_dist(lon1,lat1,lon2,lat2)

  radius=6371;
  lat1=atan(0.993277*tan(lat1/180*pi))*180/pi;
  lat2=atan(0.993277*tan(lat2/180*pi))*180/pi;
  
  temp=sin((90-lat1)/180*pi)*cos(lon1/180*pi)...
      *sin((90-lat2)/180*pi)*cos(lon2/180*pi)+sin((90-lat1)/180*pi)...
      *sin(lon1/180*pi)*sin((90-lat2)/180*pi)*sin(lon2/180*pi)...
      +cos((90-lat1)/180*pi)*cos((90-lat2)/180*pi);
  if temp > 1
      temp=1;
  end
  if temp<-1
      temp=-1;
  end
  theta=abs(acos(temp));
  dist =  theta*radius;
  
  
  
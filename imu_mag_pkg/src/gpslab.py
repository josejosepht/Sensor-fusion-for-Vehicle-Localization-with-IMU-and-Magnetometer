#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import rospy
import serial
from math import sin, pi
from std_msgs.msg import Float64
import time
import utm

from imu_mag.msg import custom

if __name__ == '__main__':
    SENSOR_NAME = "gps_sensor"
    pub= rospy.Publisher("custom_message",custom,queue_size=10)
    rospy.init_node('gps_sensor')
    serial_port = rospy.get_param('~port','/dev/ttyUSB0')
    serial_baud = rospy.get_param('~baudrate',4800)
    sampling_rate = rospy.get_param('~sampling_rate',5.0)
     
    port = serial.Serial(serial_port, serial_baud, timeout=3.)
    rospy.logdebug("Using gps sensor on port "+serial_port+" at "+str(serial_baud))
    #rospy.logdebug("Using latitude = "+str(latitude_deg)+" & atmosphere offset = "+str(offset))
    rospy.logdebug("Initializing sensor with *0100P4\\r\\n ...")
    
    sampling_count = int(round(1/(sampling_rate*0.007913)))
    rospy.sleep(0.2)        
    #line = port.readline()
     
    #latitude = latitude_deg * pi / 180.
    #depth_pub = rospy.Publisher(SENSOR_NAME+'/depth', Float64, queue_size=5)
    #pressure_pub = rospy.Publisher(SENSOR_NAME+'/pressure', Float64, queue_size=5)
    #odom_pub = rospy.Publisher(SENSOR_NAME+'/odom',Odometry, queue_size=5)
    
    rospy.logdebug("Initialization complete")
    
    rospy.loginfo("Publishing longitude and latitutde.")
        
    #odom_msg = Odometry()
    #odom_msg.header.frame_id = "odom"
    #odom_msg.child_frame_id = SENSOR_NAME
    #odom_msg.header.seq=0
    
   
    msg=custom()
    i=1
    try:
        while not rospy.is_shutdown():
            msg.header.seq=i
            line = port.readline()
            line2=line.decode('latin-1')
            #print(line2)
            if line == '':
                rospy.logwarn("DEPTH: No data")
            else:
                if line2.startswith("$GPGGA") :
                    s =line2.split(",")
                    lat = s[2]
                    lon = s[4]
                    lat_dir = s[3]
                    lon_dir = s[5]
                    utc_time = s[1]
                    alt = s[9]
 #print(lat + " lattitude" + lon + "longitude" )
 
 
#decimal degree lat=43.21583333 lon=71.76388889
#utm_data=utm.from_latlon(float(42.203177),float(-71.052450))
#print(utm_data)
#utm_data2=utm.to_latlon(328041,4689484,19,'T')
#print("\n"+str(utm_data2))
                    degrees_lat=int(float(lat)/100)
                    #print("\nDegree lat:"+str(degrees_lat))
                    minutes_lat=float(lat)-(degrees_lat*100)
                    #print("\tminutes lat:"+str(minutes_lat))
                    degrees_lon=int(float(lon)/100)
                    #print("\nDegree lon:"+str(degrees_lon))
                    minutes_lon=float(lon)-(degrees_lon*100)
                    #print("\tminutes lon:"+str(minutes_lon))
                    dd_lat= float(degrees_lat) + float(minutes_lat)/60
                    dd_lon= float(degrees_lon) + float(minutes_lon)/60 
                    if lon_dir == 'W':
                        dd_lon *= -1
                    if lat_dir == 'S':
                        dd_lat *= -1
                    print("\n"+str(dd_lat)+" "+str(dd_lon))
 
                    utm_data3=utm.from_latlon(dd_lat,dd_lon)
                    print(utm_data3)
                    msg.header.stamp=rospy.get_rostime()
                    msg.header.frame_id="GPS_Data"
                    msg.latitude=dd_lat
                    msg.longitude=dd_lon
                    msg.altitude=float(alt)
                    msg.utm_easting=utm_data3[0]
                    msg.utm_northing=utm_data3[1]
                    msg.zone=float(utm_data3[2])
                    msg.letter_field=utm_data3[3]
                    # msg.header.frame_id="GPS_Data"
                    # msg.latitude=float(1)
                    # msg.longitude=float(1)
                    # msg.altitude=float(1)
                    # msg.utm_easting=float(1)
                    # msg.utm_northing=float(1)
                    # msg.zone=float(1)
                    # msg.letter_field="T"
                    rospy.loginfo(msg)
                    pub.publish(msg)
    except rospy.ROSInterruptException:
        port.close()
    
    except serial.serialutil.SerialException:
        rospy.loginfo("Shutting down paro_depth node...")

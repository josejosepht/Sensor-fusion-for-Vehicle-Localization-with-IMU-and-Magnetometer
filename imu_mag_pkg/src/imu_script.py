
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ast import Bytes
from mmap import PROT_WRITE
import rospy
import serial
import time
import sensor_msgs
from sensor_msgs.msg import Imu
from sensor_msgs.msg import MagneticField
from scipy.spatial.transform import Rotation
from std_msgs.msg import Header

if __name__ == '__main__':
    SENSOR_NAME = "imu_sensor"
    pub1= rospy.Publisher("imu",Imu,queue_size=10)
    pub2= rospy.Publisher("mag",MagneticField,queue_size=10)
    rospy.init_node('gps_sensor')
    serial_port = rospy.get_param('~port','/dev/ttyUSB1')
    serial_baud = rospy.get_param('~baudrate',115200)
    sampling_rate = rospy.get_param('~sampling_rate',40.0)
   
    port = serial.Serial(serial_port, serial_baud, timeout=3.)
    rospy.logdebug("Using imu sensor on port "+serial_port+" at "+str(serial_baud))
    rospy.logdebug("Initializing sensor with *0100P4\\r\\n ...")
    sampling_count = int(round(1/(sampling_rate*0.007913)))
    rospy.sleep(0.2)   
    port.write(b'$VNWRG,07,40*XX')     
    #line = port.readline()
     
    rospy.logdebug("Initialization complete")
    rospy.loginfo("Publishing IMU and Magnetometer data")
            
    h1= Header()
    #h2=Header()
    msg_imu=Imu()
    msg_mag=MagneticField()
    
    i=1

    try:
        while not rospy.is_shutdown():                                 
            line = port.readline()                              #reading the serial port data
            if line == 'no data ':                                      #if there is no data, display "no data"
                rospy.logwarn("VNYMR: No data")
            else:
                if line.startswith(b'$VNYMR') :                 #parse the data from serial port only when it starts with $VNYMR
                    print(line)
                    s =line.split(b",")                         #split the string with delimiter ","
                    yaw = s[1].decode('utf-8')                  #getting required valueS
                    pitch = s[2].decode('utf-8')
                    roll = s[3].decode('utf-8')
                    print("Yaw:"+yaw+" Pitch:"+pitch+" Roll:"+roll)
                    magx = s[4].decode('utf-8')
                    magy = s[5].decode('utf-8')
                    magz = s[6].decode('utf-8')
                    print("Magx:"+magx+" Magy:"+magy+" Magz:"+magz)
                    aclx = s[7].decode('utf-8')
                    acly = s[8].decode('utf-8')
                    aclz = s[9].decode('utf-8')
                    print("Aclx:"+aclx+" Acly:"+acly+" Aclz:"+aclz)
                    gyrx = s[10].decode('utf-8')
                    gyry = s[11].decode('utf-8')
                    gyrz = s[12].decode('utf-8')
                    gyrz = gyrz[:-5]
                    print("Gyrx:"+gyrx+" Gyry:"+gyry+" Gyrz:"+gyrz)
                    rot = Rotation.from_euler('xyz', [float(roll), float(pitch), float(yaw)], degrees=True)
                    quat=rot.as_quat()
                    x=quat[0]
                    y=quat[1]
                    z=quat[2]
                    w=quat[3]
                    print("x:"+str(x)+" y:"+str(y)+" z:"+str(z)+" w:"+str(w))

                    #publish data to the message
                    h1.seq=i
                    
                    h1.stamp=rospy.get_rostime()
                    h1.stamp=rospy.get_rostime()
                    h1.frame_id="IMU DATA"
                    h1.frame_id="MAG DATA"
                    msg_imu.header=h1
                    msg_mag.header=h1
                    msg_imu.orientation.x=float(x)
                    msg_imu.orientation.y=float(y)
                    msg_imu.orientation.z=float(z)
                    msg_imu.orientation.w=float(w)
                    msg_imu.angular_velocity.x=float(gyrx)
                    msg_imu.angular_velocity.y=float(gyry)
                    msg_imu.angular_velocity.z=float(gyrz)
                    msg_imu.linear_acceleration.x=float(aclx)
                    msg_imu.linear_acceleration.y=float(acly)
                    msg_imu.linear_acceleration.z=float(aclz)
                    msg_mag.magnetic_field.x=float(magx)
                    msg_mag.magnetic_field.y=float(magy)
                    msg_mag.magnetic_field.z=float(magz)
                    pub1.publish(msg_imu)
                    pub2.publish(msg_mag)
                    i=i+1                                       #incrementing counter value for message header sequence
    except rospy.ROSInterruptException:
        port.close()
    




%Imported data from C2:C60600,E2:H60600,R2:T60600,AD2:AF60600 from imu file
%Imported data from C2:C60600,E2:G60600 from mag file
%Imported data from C2:C1518,H2:I1518 from gps file
close all
%Run only once
bag5=rosbag("imu_mag_complete.bag");
bag1=rosbag("imu_mag_stat1.bag");
bag2=rosbag("imu_mag_circles.bag");
bag3=rosbag("imu_mag_driving.bag");
bag4=rosbag("imu_mag_stat2.bag");
bsel5_2=select(bag5,"Topic","mag");
msgStruct5_2=readMessages(bsel5_2,'DataFormat','struct');
bsel1_1=select(bag1,"Topic","imu");
bsel1_2=select(bag1,"Topic","mag");
bsel1_3=select(bag1,"Topic","custom_message");
msgStruct1_1=readMessages(bsel1_1,'DataFormat','struct');
msgStruct1_2=readMessages(bsel1_2,'DataFormat','struct');
msgStruct1_3=readMessages(bsel1_3,'DataFormat','struct');
bsel2_1=select(bag2,"Topic","imu");
bsel2_2=select(bag2,"Topic","mag");
bsel2_3=select(bag2,"Topic","custom_message");
msgStruct2_1=readMessages(bsel2_1,'DataFormat','struct');
msgStruct2_2=readMessages(bsel2_2,'DataFormat','struct');
msgStruct2_3=readMessages(bsel2_3,'DataFormat','struct');
bsel3_1=select(bag3,"Topic","imu");
bsel3_2=select(bag3,"Topic","mag");
bsel3_3=select(bag3,"Topic","custom_message");
msgStruct3_1=readMessages(bsel3_1,'DataFormat','struct');
msgStruct3_2=readMessages(bsel3_2,'DataFormat','struct');
msgStruct3_3=readMessages(bsel3_3,'DataFormat','struct');
bsel4_1=select(bag4,"Topic","imu");
bsel4_2=select(bag4,"Topic","mag");
bsel4_3=select(bag4,"Topic","custom_message");
msgStruct4_1=readMessages(bsel4_1,'DataFormat','struct');
msgStruct4_2=readMessages(bsel4_2,'DataFormat','struct');
msgStruct4_3=readMessages(bsel4_3,'DataFormat','struct');
    
or_x1=cellfun(@(m) double(m.Orientation.X),msgStruct2_1);
or_y1=cellfun(@(m) double(m.Orientation.Y),msgStruct2_1);
or_z1=cellfun(@(m) double(m.Orientation.Z),msgStruct2_1);
or_w1=cellfun(@(m) double(m.Orientation.W),msgStruct2_1);
angv_x1=cellfun(@(m) double(m.AngularVelocity.X),msgStruct2_1);
angv_y1=cellfun(@(m) double(m.AngularVelocity.Y),msgStruct2_1);
angv_z1=cellfun(@(m) double(m.AngularVelocity.Z),msgStruct2_1);
linacc_x1=cellfun(@(m) double(m.LinearAcceleration.X),msgStruct2_1);
linacc_y1=cellfun(@(m) double(m.LinearAcceleration.Y),msgStruct2_1);
linacc_z1=cellfun(@(m) double(m.LinearAcceleration.Z),msgStruct2_1);
mag_x1=cellfun(@(m) double(m.MagneticField_.X),msgStruct2_2);
mag_y1=cellfun(@(m) double(m.MagneticField_.Y),msgStruct2_2);
mag_z1=cellfun(@(m) double(m.MagneticField_.Z),msgStruct2_2);
utm_easting1=cellfun(@(m) double(m.UtmEasting),msgStruct2_3);
utm_northing1=cellfun(@(m) double(m.UtmNorthing),msgStruct2_3);

or_x2=cellfun(@(m) double(m.Orientation.X),msgStruct3_1);
or_y2=cellfun(@(m) double(m.Orientation.Y),msgStruct3_1);
or_z2=cellfun(@(m) double(m.Orientation.Z),msgStruct3_1);
or_w2=cellfun(@(m) double(m.Orientation.W),msgStruct3_1);
angv_x2=cellfun(@(m) double(m.AngularVelocity.X),msgStruct3_1);
angv_y2=cellfun(@(m) double(m.AngularVelocity.Y),msgStruct3_1);
angv_z2=cellfun(@(m) double(m.AngularVelocity.Z),msgStruct3_1);
linacc_x2=cellfun(@(m) double(m.LinearAcceleration.X),msgStruct3_1);
linacc_y2=cellfun(@(m) double(m.LinearAcceleration.Y),msgStruct3_1);
linacc_z2=cellfun(@(m) double(m.LinearAcceleration.Z),msgStruct3_1);
mag_x2=cellfun(@(m) double(m.MagneticField_.X),msgStruct3_2);
mag_y2=cellfun(@(m) double(m.MagneticField_.Y),msgStruct3_2);
mag_z2=cellfun(@(m) double(m.MagneticField_.Z),msgStruct3_2);
utm_easting2=cellfun(@(m) double(m.UtmEasting),msgStruct3_3);
utm_northing2=cellfun(@(m) double(m.UtmNorthing),msgStruct3_3);

a1_1=cellfun(@(m) double(m.Header.Stamp.Sec),msgStruct2_1);
a1_2=cellfun(@(m) double(m.Header.Stamp.Sec),msgStruct2_2);
a1_3=cellfun(@(m) double(m.Header.Stamp.Sec),msgStruct2_3);
x=cellfun(@(m) double(m.Header.Stamp.Sec),msgStruct3_1);
y=cellfun(@(m) double(m.Header.Stamp.Sec),msgStruct3_2);
z=cellfun(@(m) double(m.Header.Stamp.Sec),msgStruct3_3);
b1_1=[a1_1;x];
b1_2=[a1_2;y];
b1_3=[a1_3;z];
imu_time =b1_1;
mag_time = b1_2;
gps_time = b1_3;
imu_orientation = [or_x1(:,:) or_y1(:,:) or_z1(:,:) or_w1(:,:);or_x2(:,:) or_y2(:,:) or_z2(:,:) or_w2(:,:)];
imu_ang_vel = [mag_x1 mag_y1 mag_z1;mag_x2 mag_y2 mag_z2];
imu_lin_acc = [linacc_x1 linacc_y1 linacc_z1;linacc_x2 linacc_y2 linacc_z2];

mag_measurements =[mag_x1 mag_y1 mag_z1;mag_x2 mag_y2 mag_z2];

gps_measurements = [utm_easting1 utm_northing1;utm_easting2 utm_northing2];

sz = 40;

%Bias subtraction from acceleration
bias = [mean(imu_lin_acc(:,1)) mean(imu_lin_acc(:,2)) mean(imu_lin_acc(:,3))];
imu_lin_acc_unbiased = [imu_lin_acc(:,1)-bias(1) imu_lin_acc(:,2)-bias(2) imu_lin_acc(:,3)-bias(3)];

gps_time_lag = imu_time(1)-gps_time(1);
gps_time_matched = gps_time+gps_time_lag;


imu_time_real = imu_time;
imu_orientation_real = imu_orientation;
imu_ang_vel_real = imu_ang_vel;
imu_lin_acc_real = imu_lin_acc_unbiased;
imu_lin_acc_real = highpass(imu_lin_acc_real,15,40);

mag_time_real = mag_time;
mag_measurements_real = mag_measurements;
%mag_measurements_real=lowpass(mag_measurements_real,(1/40));

gps_time_real = gps_time_matched;
gps_measurements_real = gps_measurements;
mag_y5 = cellfun(@(m) double(m.MagneticField_.Y),msgStruct5_2);
mag_x5 = cellfun(@(m) double(m.MagneticField_.X),msgStruct5_2);
mag_z5 = cellfun(@(m) double(m.MagneticField_.Z),msgStruct5_2);
figure
plot([utm_easting1;utm_easting2],[utm_northing1;utm_northing2])
title("Plot of GPS Data")
xlabel("Utm Easting")
ylabel("Utm Northing")

figure
subplot(2,2,1)
plot([mag_x1;mag_x2],[mag_y1;mag_y2])
title("Magnetometer readings(X and Y)")
xlabel("Magnetometer X reading(Gauss)")
ylabel("Magnetometer Y reading(Gauss)")
subplot(2,2,2)
plot([mag_x1;mag_x2])
title("Magnetometer X reading w.r.t time")
xlabel("Magnetometer X reading(Gauss)")
ylabel("Time")
subplot(2,2,3)
plot([mag_y1;mag_y2])
title("Magnetometer Y reading w.r.t time")
xlabel("Magnetometer Y reading(Gauss)")
ylabel("Time")
subplot(2,2,4)
plot([mag_z1;mag_z2])
title("Magnetometer Z reading w.r.t time")
xlabel("Magnetometer Z reading(Gauss)")
ylabel("Time")

figure
plot(mag_x1,mag_y1)
title("Magnetometer readings plot X vs Y for circular motion data")
xlabel("Magnetometer X reading(Gauss)")
ylabel("Magnetometer Y reading(Gauss)")
grid on

%Hard Iron effects
a=mean(mag_x1);
b=mean(mag_y1);
figure
mag_x1_hcorrected=mag_x1-a;
mag_y1_hcorrected=mag_y1-b;
mag_x2_hcorrected=mag_x2-a;
mag_y2_hcorrected=mag_y2-b;
plot(mag_x1_hcorrected,mag_y1_hcorrected)
title("Magnetometer readings plot X vs Y after Hard iron effect correction for circular motion data")
xlabel("Magnetometer X reading(Gauss)")
ylabel("Magnetometer Y reading(Gauss)")
grid on

%Soft iron effect
r = ((mag_x1_hcorrected.^2) + (mag_y1_hcorrected.^2)).^0.5;
[m_v,m_I] = min(r);
[r_v,r_I] = max(r);

theta = asind(mag_y1_hcorrected(r_I)/r_v);
R = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
A=[mag_x1_hcorrected mag_y1_hcorrected;mag_x2_hcorrected mag_y2_hcorrected];
mag_corrected_values_real = R*A';
temp_mag_circle_corrected=[mag_x1_hcorrected mag_y1_hcorrected];
v_1 = R*temp_mag_circle_corrected';

scale_factor = m_v/r_v;
mag_corrected_values_real(1,:) = mag_corrected_values_real(1,:).*scale_factor;
v_1(1,:) = v_1(1,:).*scale_factor;
theta = -theta;
R = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
mag_corrected_values_real = R*mag_corrected_values_real;
v_1 = R*v_1;

figure 
axis equal;
v_1=lowpass(v_1,0.9);
plot(v_1(1,:),v_1(2,:));
title('Plot of magnetometer readings for circles after adjustments for Soft Iron effects'); 
xlabel('X-readings(Gauss)'); 
ylabel('Y-readings(Gauss)');
grid on

figure
axis equal; axis square;
plot(mag_corrected_values_real(1,:),mag_corrected_values_real(2,:));
title('Plot of magnetometer readings after adjustments for Soft and Hard Iron effects'); xlabel('X-readings'); ylabel('Y-readings');

mag_corrected_values_real = mag_corrected_values_real';
yaw_mag = atand(mag_corrected_values_real(:,1)./mag_corrected_values_real(:,2));

imu_time_diff = (imu_time_real-imu_time_real(1))./10^9;
yaw_ang_vel_rad = cumtrapz([angv_z1;angv_z2]);
yaw_ang_vel_deg = wrapTo180(rad2deg(yaw_ang_vel_rad));
%yaw_ang_vel_deg_int = cumtrapz(yaw_ang_vel_deg);
eul = quat2eul([[or_w1;or_w2] [or_x1;or_x2] [or_y1;or_y2] [or_z1;or_z2]]);
yaw_eul = wrapTo180(rad2deg(eul(:,1)));

figure
hold on
plot(unwrap(yaw_mag)./7);
plot(unwrap(-yaw_ang_vel_deg)./20);
plot(unwrap(yaw_eul)./7);
title('Plot of yaw calculated from magnetometer, gyro, and orientation'); 
legend({'Magnetometer','Gyro','Orientation'});
xlabel('time'); 
ylabel('yaw angle (degrees)');
hold off

a = 0.4;
comp_fil_yaw = (yaw_mag.*a)+(yaw_ang_vel_deg(1:39017,1).*(1-a));

figure
plot(unwrap(comp_fil_yaw));
hold on
plot(unwrap(yaw_eul));
title('Plot of yaw calculated using Complimentary filter and from orientation'); 
legend({'Complimentary Filter','Orientation'}); 
xlabel('time'); 
ylabel('yaw angle (degrees)');
hold off


%Part 3.2:
imu_velocity = cumtrapz(imu_time_diff, imu_lin_acc_real(:,1));
imu_velocity=imu_velocity.*10^10;
% gps_vel_bias = [utm_easting2(1) utm_northing2(1)]
% gps_velocity = zeros(size(gps_vel_bias,1));
% for c = 1:size(gps_vel_bias,1)-1
%      temp = gps_vel_bias(c+1,:)-gps_vel_bias(c,:);
%      gps_velocity(c) = ((temp(1)^2 +temp(2)^2)^0.5)/((gps_time(c+1)-gps_time(c)));
% end

gps_velocity = (zeros(size(gps_measurements_real,1)));
%gps_velocity = gps_velocity -gps_vel_bias;
for c = 1:size(gps_measurements_real,1)-1
     temp = gps_measurements_real(c+1,:)-gps_measurements_real(c,:);
     gps_velocity(c) = ((temp(1)^2 +temp(2)^2)^0.5)/((gps_time(c+1)-gps_time(c)));
end

figure
plot(gps_time_real,gps_velocity);
hold on
plot(imu_time_real(1:size(imu_velocity)), imu_velocity);
title('Plot of velocities from GPS and from IMU'); legend('','IMU Velocity');xlabel('time'); ylabel('velocity');
hold off

y_dot_dot_calculated = imu_ang_vel_real(:,3).*imu_velocity(1:39017,:);

figure
plot(imu_time_real, imu_lin_acc_real(:,2));
hold on
plot(imu_time_real(1:39017), y_dot_dot_calculated);
title('Plot of acceleration of y_dot_dot_observed and w*X_dot'); legend('w*X_dot','y_dot_dot_observed');xlabel('time'); ylabel('acceleration');
hold off



v_e = imu_velocity.*cos(deg2rad(yaw_eul));
v_n = imu_velocity.*sin(deg2rad(yaw_eul));

x_n = cumtrapz(mag_time_real./10^9,v_n(1:39017));
x_e = cumtrapz(mag_time_real./10^9,v_e(1:39017));

figure
plot(gps_time_real,gps_measurements_real(:,1)-gps_measurements_real(1,1));
hold on
plot(imu_time_real(1:39017),x_e*(-(10^9.3)));
title('Plot of Displacements obtained from GPS and IMU (X-axis)'); legend('GPS','IMU');xlabel('time'); ylabel('X-axis/UTM Easting');
hold off

figure
plot(gps_time_real,gps_measurements_real(:,2)-gps_measurements_real(1,2));
hold on
plot(imu_time_real(1:39017),x_n*((10^9)));
title('Plot of Displacements obtained from GPS and IMU (Y-axis)'); legend('GPS','IMU');xlabel('time'); ylabel('Y-axis/UTM Northing');
hold off
%{
th = 30;
R = [cosd(th) sind(th); -sind(th) cosd(th)];
scale = 1.5;

x_e = -1.*x_e;
% x_n = -1.*x_n;

imu_trajectory = scale.*(R*[x_e x_n]');

figure
plot(imu_trajectory(1,1:39017).*((1/1.75)*(10^9)),imu_trajectory(2,1:39017).*((7/0.5)*10^8),'b');
hold on
plot(gps_measurements_real(:,1)-gps_measurements_real(1,1),gps_measurements_real(:,2)-gps_measurements_real(1,2),'r');
title('Plot of trajectories obtained from GPS and IMU'); legend('IMU','GPS');xlabel('X-axis/UTM Easting'); ylabel('Y-axis/UTM Northing');
%}

%th = 200;
th = 25;
R = [cosd(th) sind(th); -sind(th) cosd(th)];
scale = 0.9;

x_e = -1.*x_e;
% x_n = -1.*x_n;

imu_trajectory = scale.*(R*[x_e x_n]');

figure
plot(imu_trajectory(1,1:39017).*((1/1.75)*(10^9)),imu_trajectory(2,1:39017).*((7/0.5)*10^8),'b');
hold on
plot(gps_measurements_real(:,1)-gps_measurements_real(1,1),gps_measurements_real(:,2)-gps_measurements_real(1,2),'r');
title('Plot of trajectories obtained from GPS and IMU'); legend('IMU','GPS');xlabel('X-axis/UTM Easting'); ylabel('Y-axis/UTM Northing');

x=(imu_lin_acc_real(1:39017,2)-y_dot_dot_calculated);
ang_acc=imu_velocity(2:39018,1)-imu_velocity(1:39017,1);
ang_acc=ang_acc*40;
for i=1:39017
    if ang_acc(i)~=0
        x_c(i)=x(i)/ang_acc(i);
    end
end
mean(x_c)

%scratch code for debugging plots
% th = -90;
% R = [cosd(th) sind(th); -sind(th) cosd(th)];
% scale =1.5;
% 
% x_e = imu_trajectory(1,22571:39017);
% x_n = imu_trajectory(2,22571:39017);
% 
% imu_trajectory = scale.*(R*[x_e' x_n']');
% plot((imu_trajectory(1,:).*((1/1.75)*(10^9)))+860.325,(imu_trajectory(2,:).*((7/0.5)*10^8))-548.296,'b');
% 
% th = -90;
% R = [cosd(th) sind(th); -sind(th) cosd(th)];
% scale =1.5;
% 
% x_e = imu_trajectory(1,7010:16447);
% x_n = imu_trajectory(2,7010:16447);
% 
% imu_trajectory = scale.*(R*[x_e' x_n']');
% plot((imu_trajectory(1,7010:9438).*((1/1.75)*(10^9)))+860.325,(imu_trajectory(2,7010:9438).*((7/0.5)*10^8))-548.296,'b');
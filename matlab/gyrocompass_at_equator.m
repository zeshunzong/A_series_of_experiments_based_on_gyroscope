close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS CODE AIMS AS SIMULATE HOW A GYROCOMPASS WORKS
% FOR MATH-UA 395 PROJECT 1
% AUTHOR: ZESHUN ZONG AND YIFENG CHEN
% CURRENTLY NO FRICTION IS ADDED IN THIS MODEL, SO ONLY OSCILLATION CAN
% BE OBSEREVED. IN FURTHER STUDIES, THE PERIOD OF OSCILLATION CAN BE
% STUDIED. ALSO, FRICTION CAN BE ADDED TO OBTAINED A STABLE GYROCOMPASS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE CONSTRUCTION OF A GYROCOMPOSS ON THE EARTH CAN BE DECOMPOSED INTO
% THREE PARTS
% 1) BUILD AN EARTH
% 2) BUILD A TRACK AT SOME PARTICULAR POINT ON THE EARTH
% 3) ASSEMBLE A GYROSCOPE IN THIS TRACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) BUILD AN EARTH
% parameters to build an earth:
NDIM = 3;
Rearth = 100; % radius of earth
nseg = 12; % number of segments on a latitude line
ncir = 7; % number of latitde lines, one equator, (ncir-1)/2 lines on each hemisphere

% we first need to build an earth
% we do so by first using spherical coordinates (R, theta, phi)
% we divide each latitude line into nseg segments
thetas = 0 : 2*pi/nseg : 2*pi - 0.001;
% we cut the ball horiontally to get ncir circles
phis = pi/(ncir+1) : pi/(ncir+1) : pi-0.001;
% so in total ncir circles, each circle has nseg points on it. ncir * nseg
% nodes all together.
% the nodes are indexed from the bottom to the top, one circle after the
% other

Xearth = zeros(ncir*nseg, NDIM); % positions of nodes that represent the earth
for i = 1:ncir
    % loop through the circles
    for j = 1:nseg
        % loop through the points on the circles
        % x = r*cos(theta)*sin(phi); y = r*sin(theta)*sin(phi); z =
        % r*cos(phi).
        Xearth(nseg*(i-1)+j,:) = [Rearth*cos(thetas(j))*sin(phis(i)), ...
            Rearth*sin(thetas(j))*sin(phis(i)), Rearth*cos(phis(i))];
    end
end


% how many links are there? on each circle there are nseg links, 
% so nseg*ncir links here.
% jje and kke are jj and kk for the earth, denoting the staring point and
% ending point for each link
jje = zeros(ncir*nseg,1);
kke = zeros(ncir*nseg,1);
% draw latitudes
% this should fill in the first nseg*ncir links
for i = 1:ncir
    jje((i-1)*nseg+1 : i*nseg) = (i-1)*nseg+1 : i*nseg;
    kke((i-1)*nseg+1 : i*nseg) = [(i-1)*nseg+2 : i*nseg, (i-1)*nseg+1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) PICK A POINT ON THE EQUATOR, INSTALL A CYCLIC TRACK HERE 
% parameters to build the track
Rbase = 20; % the radius of the track

% HERE WE SPECIFY THAT THE POINT IS THE FIRST POINT ON THE EQUATOR.
% PARTICULARLY, THE POSITION VECTOR OF THIS POINT W.R.T. THE CENTER OF
% EARTH IS [Rearth, 0, 0]
gyropos = Xearth(((ncir-1)/2 * nseg)+1,:);

% the vector "gyropos" points from the origin (center of earth) to the
% gyrocompass (center of gyroscope). This vector, which is also radius, is 
% perpendicular to the plan tangent to the earth at gyropos

% WE MANUALLY SPECIFY TWO UNIT VECTORS TO GENERATE THE TRACK
% unitvec1 WOULD ALSO BE THE INITIAL DIRECTION OF THE AXIS OF GYROSCOPE,
% HERE WE SET ITS INITIAL DIRECTION TO BE PARALLEL WITH THE EQUATOR

%unitvec1 = [0 1 0];
%unitvec2 = [0 0 1];

unitvec1 = [0 sqrt(2)/2 sqrt(2)/2];
unitvec2 = [0 sqrt(2)/2 -sqrt(2)/2];

% suppose we want 61 points on the tangent base. 60 points uniformly
% distributed on the circumference, plus one more point at the center of
% track (i.e. gyropos)

thetas = 0: 2*pi/60: 2*pi-0.01; % angles to pin down positions of points

Xbase = zeros(61, 3); % positions of points on the track
for i = 1:60
    % the 60 points on the circumference are given by linear combinations
    % of the two unit vectors 
    Xbase(i,:) = gyropos + Rbase .* cos(thetas(i))*unitvec1 + ...
        Rbase .* sin(thetas(i))*unitvec2;
end
% the last point is the center of the track
Xbase(61,:) = gyropos;

% there are also 60 links, and 2 more links showing unitvec1 and unitvec2
%jjbase = [1:60, 61, 61];
%kkbase = [2:60, 1, 1, 46];

jjbase = [1:60];
kkbase = [2:60, 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) ASSEMBLE A GYROSCOPE THAT FITS IN THE TRACK

% let us attach a gyroscope to this track, we do so by first build a
% gyroscope centered at the origin, then change its orientation and shift
% its center of mass
% this gyroscope should have a radius of Rbase

% number of nodes on the circumference of the wheel
n = 30;

% total of (n+3) points
% n points on the circle, three points on the axis. two endpoints, the
% other is the center of the wheel
% the n points on the circle are indexed from 1 to n
% center: n+1
% top: n+2
% bottom: n+3

% total of 2n+2 links
% n spokes, indexed from 1 to n, connect from the center to the points on
% the wheel
% n edges between the n points on the wheel, indexed from n+1 to 2n
% center to the top: link (2n+1)
% center to the bottom: link (2n+2)
% initialize arrays
num_nodes = n+3;
num_links = 2*n+2;
jj = zeros(num_links,1);
kk = zeros(num_links,1);
X = zeros(num_nodes, 3);
U = zeros(num_nodes, 3);

% STEP ONE: build vertical structure
center = [0, 0, 0];
for k = 1:n
    theta = 2*pi*k/n;
    X(k,:) = center + Rbase*[cos(theta), sin(theta), 0];
end
X(n+1,:) = center;
% set position of the two other endpoints
X(n+2,:) = center + [0,0,Rbase];
X(n+3,:) = center - [0,0,Rbase]; % bottom should be at origin

% STEP TWO: rotate the structure
% we want to rotate the gyroscope so that its rotating axis is parallel to
% the direction of unitvec1.
% we do so by (change of coordinates) -> (rotation) -> (change back)
direc = Xbase(14,:) - Xbase(61,:);
avec = [direc(1), direc(2), 0]; % projection of unitvec1 on xOy plane
avec = avec ./ norm(avec); % normalize it, this is the new x direction
cvec = [0 0 1]; % we still use the original z direction
bvec = cross(cvec, avec); % construct the other direction by cross product
% rotate about the b axis
phi = acos(unitvec1(3));
rotateMat = [cos(phi), 0, sin(phi); 0, 1, 0; -sin(phi), 0, cos(phi)];
changeOfCoordMat = [avec; bvec; cvec];

% the action is P^-1 R P
action = changeOfCoordMat \ (rotateMat * changeOfCoordMat);

% STEP THREE, shift the gyroscope
% center should be at gyropos, instead of [0 0 0]

for k = 1:n+3
    X(k,:) = (action * X(k,:)')';
    X(k,:) = X(k,:) + gyropos;
end

% naming index sets for the links
spokes = 1:n;
rimlinks = (n+1):2*n;
center2top = 2*n+1;
center2bottom = 2*n+2;

% build structure by creating links
jj(spokes) = 1:n;
kk(spokes) = n+1;
jj(rimlinks) = 1:n;
kk(rimlinks) = [2:n,1];
jj(center2top) = n+1;
jj(center2bottom) = n+1;
kk(center2top) = n+2;
kk(center2bottom) = n+3;

% mass of the nodes on the gyroscope
M = ones(num_nodes,1).*0.0001;

% initialize the position of the center of mass
Xcm = (sum((M.*X))./sum(M))';
% initialize velocity center of mass
Ucm = (sum((M.*U))./sum(M))';

% initialize Xtwiddle, which is the position vectors relative to the center
% of the mass
Xtwiddle = X - Xcm';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS TO BE CHANGED


% the initial angular momentum of gyroscope
L = action * [0; 0; 10000];
% it is also rotated so that it is parallel with the rotating axis


dt = 1e-8;
%end_time = 0.02 * 4e-1; 
end_time = 0.01 * 4e-1; 
timevec = 0:dt:end_time;

% the angular velocity of the rotation of earth
earth_omega = 1000;

% magnitude of the correcting force    
stiff = 50000;


% add friction to stabilize the axis of gyroscope. The friction is induced
% by the gliding of gyroscope's axis in the track
%friction_level = 15000;
friction_level = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Draw the initial structure

%{
figure(1);
% plot the earth
x = [Xearth(jje,1) Xearth(kke,1)];
y = [Xearth(jje,2) Xearth(kke,2)];
z = [Xearth(jje,3) Xearth(kke,3)];
plot3(x',y',z','linewidth',1)
axis equal
hold
% plot the radius at this point, this is the normal vector
plot3([gyropos(1) , 0], [gyropos(2) , 0], [gyropos(3) , 0],  'linewidth', 1)
% plot the track
x = [Xbase(jjbase,1) Xbase(kkbase,1)];
y = [Xbase(jjbase,2) Xbase(kkbase,2)];
z = [Xbase(jjbase,3) Xbase(kkbase,3)];
plot3(x',y',z','linewidth',2)


% plot the gyroscope
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
z = [X(jj,3) X(kk,3)];
plot3(x',y',z','linewidth',3)
view(50,7)



%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIMULATION BEGINS HERE

% the angle that the earth rotates in a dt period
theta_per_dt = earth_omega * dt;
% the rotation of earth is governed by a rotation matrix about z axis
rotate_earth = [cos(theta_per_dt), -sin(theta_per_dt), 0; ...
    sin(theta_per_dt), cos(theta_per_dt), 0; 0, 0, 1];

% initialize some variables that will be used in the loop
P = eye(3) - cross(unitvec1, unitvec2) * (cross(unitvec1, unitvec2)');
Xtoplast = X(n+2,:);
Xbotlast = X(n+3,:);


% initialize stuff to make a movie
v = VideoWriter('gyro_EQ_init_NE');
open(v);

figure(2);

Angles = zeros(floor(length(timevec)/700)+1, 1);

for t = 1:length(timevec)

    % this just displays the current simulation time in the matlab window
    % disp(timevec(t));	
    
    
    % we view the track as fixed on the earth.
    % so they are rotated by the rotation matrix simultaneously 
    for i = 1:length(Xearth)
        Xearth(i,:) = (rotate_earth * (Xearth(i,:)'))';
    end
    for i = 1:length(Xbase)
        Xbase(i,:) = (rotate_earth * (Xbase(i,:)'))';
    end
    % also rotate the two unit vectors that span the tangent plane
    unitvec1 = (rotate_earth * (unitvec1'))';
    unitvec2 = (rotate_earth * (unitvec2'))';
    % the projection matrix that projects the axis onto the track is given
    % by P = (I - (u cross v)*(u cross v)^T)
    P = eye(3) - cross(unitvec1, unitvec2) * (cross(unitvec1, unitvec2)');
    
    % the current projected position on the track, minus the previous
    % projected position on the track, is the approximated gliding along
    % the track
    %Xtopdisplacement = (P * (X(n+2,:)'))' - Xtoplast;
    %Xbotdisplacement = (P * (X(n+3,:)'))' - Xbotlast;
    
    Xtopdisplacement = X(n+2,:) - Xtoplast;
    Xbotdisplacement = X(n+3,:) - Xbotlast;
    
    
    
    % here we assume the center of gyroscope is fixed, always centered at
    % the track. 
    trans_vec = Xbase(61,:) - X(n+1,:);
    Xcm = Xcm + trans_vec';


    % compute moment of inertia tensor
    I = zeros(NDIM, NDIM);
    for k = 1:num_nodes
        I = I + M(k).*((norm(Xtwiddle(k,:))^2).*eye(NDIM) - ...
            Xtwiddle(k,:)'*Xtwiddle(k,:) );
    end  
    
    % compute the angular velocity
    Omega = I\L;
    
    % normalize angular velocity if it is nonzero
    if(norm(Omega) > 100*eps)
         unit_Omega = Omega/norm(Omega);
         Omega_cross = [0 -Omega(3) Omega(2); Omega(3) 0 -Omega(1);...
             -Omega(2) Omega(1) 0];
         P_Omega = unit_Omega*unit_Omega';
         Xtwiddle = (P_Omega*(Xtwiddle') + ...
             cos(norm(Omega)*dt).*(eye(NDIM) - P_Omega)*(Xtwiddle') + ...
             sin(norm(Omega)*dt).*(Omega_cross*(Xtwiddle'))./norm(Omega) )';
    end

    % compute net force and net torque   
    net_force = zeros(NDIM,1);
    net_torque = zeros(NDIM,1);
    
    % Due to the earth's rotation, it will have some elevation
    % This will be corrected by forces from the track, by basically pulling
    % the two endpoints back to the track 
    [Ftop, Fbottom] = force_by_track(Xbase, X, n, stiff);
    
    
    % add friction here, friction is proportional to the velocity, pointing
    % into the opposite direction. Here we approximate the velocity by
    % difference of position vector X(t+1) - X(t), so the friction is 
    % -k * (X(t+1) - X(t)). Friction acts on both end points
    
    Ftop = Ftop - friction_level * Xtopdisplacement;
    Fbottom = Fbottom - friction_level * Xbotdisplacement;
    
    external_force = zeros(num_nodes, NDIM);
    external_force(n+2,:) = Ftop;
    external_force(n+3,:) = Fbottom;

    
    for l = 1:num_nodes
    	net_force = net_force + external_force(l,:)';
        net_torque = net_torque + cross(Xtwiddle(l,:)', external_force(l,:)');
    end
    
 
    % update the position and velocity for the center of mass
    Ucm = Ucm + (dt/sum(M)).*net_force;
    Xcm = Xcm + dt.*Ucm;

    % update the angular momentum
    L = L + dt.*net_torque;
    

    % update positions of individual masses
    X = Xtwiddle + Xcm';
    
    
    % here we track the (corresponding) positions of the two endpoints on
    % the track so that we can use later to add friction
 
    Xtoplast = (P * X(n+2,:)')';
    Xbotlast = (P * X(n+3,:)')';

    
    mark = mod(t,700);
    
    if mark < 1
        
        
        % plot the earth
        xe = [Xearth(jje,1) Xearth(kke,1)];
        ye = [Xearth(jje,2) Xearth(kke,2)];
        ze = [Xearth(jje,3) Xearth(kke,3)];
        plot3(xe',ye',ze','linewidth',1)
        axis equal
        hold on;
        % plot the track
        xb = [Xbase(jjbase,1) Xbase(kkbase,1)];
        yb = [Xbase(jjbase,2) Xbase(kkbase,2)];
        zb = [Xbase(jjbase,3) Xbase(kkbase,3)];
        plot3(xb',yb',zb','linewidth',2)
        % plot the gyroscope
        xg = [X(jj,1) X(kk,1)];
        yg = [X(jj,2) X(kk,2)];
        zg = [X(jj,3) X(kk,3)];
        plot3(xg',yg',zg','linewidth',3)
        %az = 250;
        az = 95 + 55 * t * theta_per_dt;
        el = 10;
        view(az, el);
        hold off;
        
        frame = getframe(gcf);
        writeVideo(v,frame);
        
        
        % we also want to record the deviating angle each time we plot the
        % gyrocompass. Suppose the vector pointing north is n, the axis is
        % a, then cos(angle) = dot(a, n)/ (norm(a)*norm(n))
        % here notice that a should be the axis projected onto the plane
        
        %axis_direction = (P * (X(n+2,:) - X(n+1,:))')';
        axis_direction = X(n+2,:) - X(n+1,:);
        cosangle = dot([0 0 1], axis_direction)/(norm([0 0 1])*norm(axis_direction));

        Angles(floor(t/700),:) = acos(cosangle);
        
    end

    %axis equal
    %view([0,5])
    %xlim([-3 3])
    %ylim([-1 1])
    %zlim([-2 4])
    pause(0.0000001)
 
  
    % stop simulation if it blows up.
    if(norm(X)> 1e10)
      break
    end

end

close(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SIMULATION

function [Ftop, Fbottom] = force_by_track(Xbase, X, n, stiff)


    % the tangent plane is spanned by two unit vectors on the plane
    vec1 = Xbase(1,:) - Xbase(61,:);
    vec2 = Xbase(46,:) - Xbase(61,:);
    normal_vec = cross(vec1, vec2);
    % the equation of the plane is given by 
    % n1(x-a) + n2(y-b) + n3(z-c) = 0, i.e. 
    % n1*x + n2*y + n3*z - a*n1 - b*n2 - c*n3 = 0
    % n1*x + n2*y + n3*z + (-a*n1-b*n2-c*n3) = 0. Denote the coeff A,B,C,D
    % the distance of a point to the plane is given by 
    % (Ax+By+Cz+D)/sqrt(A^2 + B^2 + C^2)
    % if positive, on the same side with normal vector
    
    % the point that is always on the center is Xbase(61,:)
    % compute the coefficients
    A = normal_vec(1);
    B = normal_vec(2);
    C = normal_vec(3);
    D = -Xbase(61,1)*A - Xbase(61,2)*B - Xbase(61,3)*C;
    top_to_plane = (dot([A B C], X(n+2,:))+D)/sqrt(A^2 + B^2 + C^2);
    bot_to_plane = (dot([A B C], X(n+3,:))+D)/sqrt(A^2 + B^2 + C^2);
    
    % force should be propotional tp the distance

    Ftop = - top_to_plane * stiff * normal_vec;
    Fbottom = - bot_to_plane * stiff * normal_vec;

end


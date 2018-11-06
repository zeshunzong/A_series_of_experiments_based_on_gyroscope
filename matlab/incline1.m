close all; clear all; clc;


% build a gyroscope, i.e. a wheel plus an axis
% we no longer let the axis be vertically placed. Instead, we give the
% gyroscope an initial angle.

NDIM = 3;
% set numerical parameters
% we fix the total amount of time to be 0.01 second, and we will allow dt
% to change in order to check numerical stability
dt = 1e-7;
end_time = 4e-1; 
timevec = 0:dt:end_time;

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
% n spokes, indexed from 1 to n, connect from the center to the points on the
% wheel
% n edges between the n points on the wheel, indexed from n+1 to 2n
% center to the top: link (2n+1)
% center to the bottom: link (2n+2)
% initialize arrays
num_nodes = n+3;
num_links = 2*n+2;
jj = zeros(num_links,1);
kk = zeros(num_links,1);
X = zeros(num_nodes, NDIM);
U = zeros(num_nodes, NDIM);

% parameters for the physical structure of gyroscope
angle = pi/12; % default to be 15 degrees
bin_length = 0.8; % length from the center to the bottom/top
RB = 1; % radius of wheel


% initially we let the axis lie in the xz plane, deviating "angle" degrees
% from the z axis. To acheive this, we first build a vertical gyroscope,
% and then use a rotation matrix to rotate the gyroscope about the y axis
% (so that the axis of the gyroscope will fall in xOz plane)

% assume each node has mass 1
M = ones(num_nodes,1); 

% magnitude of the gravitational force
G = 1000000; 

% the rotation matrix, see https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
rotation_along_y = [cos(angle), 0, sin(angle);0,1,0;-sin(angle), 0, cos(angle)];

% set positions of nodes around the circumference of the wheel in two
% steps, first build the vertical structure, then rotate it.

% STEP ONE: build vertical structure
center = [0, 0, bin_length];
for k = 1:n
    theta = 2*pi*k/n;
    X(k,:) = center + RB*[cos(theta), sin(theta), 0];
end
X(n+1,:) = center;
% set position of the two other endpoints
X(n+2,:) = center + [0,0,bin_length];
X(n+3,:) = center - [0,0,bin_length]; % it should be that the bottom is at the origin

% STEP TWO: rotate the structure, by simply multipying the rotation matrix
for k = 1:n+3
    X(k,:) = (rotation_along_y * X(k,:)')';
end

original_bottom = X(n+3,:); % this position will be recored so that later we can fix this point by stiff springs


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

%{
% a section of code to plot the initial structure
figure(1);
x = [X(jj,1) X(kk,1)];
y = [X(jj,2) X(kk,2)];
z = [X(jj,3) X(kk,3)];
plot3(x',y',z','linewidth',4)
axis equal
view([0,6])
%}


% initialize the position of the center of mass
Xcm = (sum((M.*X))./sum(M))';
% initialize velocity center of mass
Ucm = (sum((M.*U))./sum(M))';

% initialize Xtwiddle
%Xtwiddle = zeros(num_nodes,NDIM);
Xtwiddle = X - Xcm';


% make a gravitational force
F_gravity = zeros(num_nodes, NDIM);
F_gravity(:,1) = 0;
F_gravity(:,2) = 0;
F_gravity(:,3) = -G.*M; %M is a vector denoting the mass of every node, while G is some constant, say 9.8
% F_grabity is a matrix, recording the gravitational force on every node.


% initialize the angular momentom. We let the wheel spin
L = (rotation_along_y * [0; 0; 100000]);



% here is the timestep loop
bottom_position = zeros(length(timevec), NDIM);
for t = 1:length(timevec)

    % this just displays the current simulation time in the matlab window
    % disp(timevec(t));	

    % compute moment of inertia tensor
    I = zeros(NDIM, NDIM);
    for k = 1:num_nodes
	I = I + M(k).*((norm(Xtwiddle(k,:))^2).*eye(NDIM) - Xtwiddle(k,:)'*Xtwiddle(k,:) );
    end     
    
    % compute the angular velocity
    Omega = I\L;
    
    % normalize angular velocity if it is nonzero
    if(norm(Omega) > 100*eps)
         unit_Omega = Omega/norm(Omega);
         Omega_cross = [0 -Omega(3) Omega(2); Omega(3) 0 -Omega(1); -Omega(2) Omega(1) 0];
         P_Omega = unit_Omega*unit_Omega';
         Xtwiddle = (P_Omega*(Xtwiddle') + cos(norm(Omega)*dt).*(eye(NDIM) - P_Omega)*(Xtwiddle') + sin(norm(Omega)*dt).*(Omega_cross*(Xtwiddle'))./norm(Omega) )';
    end

    % compute net force and net torque   
    net_force = zeros(NDIM,1);
    net_torque = zeros(NDIM,1);
    
    % one way to exert normal force is to make a table, the normal force is
    % pointing vertically upward, whose magnitude is proportional to how
    % deep the point is "inside" the table.
    normal_force = F_ground(X);
 
    
    %{
    % another way to exert normal force is to use springs
    % try to add three really stiff springs to fix the bottom point
    stiffness = 500000000;
    fixing_force = (original_bottom - X(n+3,:)).*stiffness;
    % we want more stiffness on the vertical axis
    fixing_force(3) = fixing_force(3)*2;
    normal_force = zeros(num_nodes, NDIM);
    normal_force(n+3,:)= fixing_force;
    %}
  

    for l = 1:num_nodes
    	net_force = net_force + F_gravity(l,:)' + normal_force(l,:)';
        net_torque = net_torque + cross(Xtwiddle(l,:)', F_gravity(l,:)' + normal_force(l,:)');
    end
 
    % update the position and velocity for the center of mass
    Ucm = Ucm + (dt/sum(M)).*net_force;
    Xcm = Xcm + dt.*Ucm;

    % update the angular momentum
    L = L + dt.*net_torque; 

    % update positions of individual masses
    X = Xtwiddle + Xcm';
    
    %X = horizontal_translate(X);

    % store the position of centroid
    bottom_position(t,:) = X(n+3,:);

    % plot the current position of the ball
    
    mark = mod(t,300);
    if mark < 1
        figure(2);
        x = [X(jj,1) X(kk,1)];
        y = [X(jj,2) X(kk,2)];
        z = [X(jj,3) X(kk,3)];
        plot3(x',y',z','linewidth',4)

        axis equal
        view([0,5])
        xlim([-3 3])
        ylim([-1 1])
        zlim([-2 4])
        pause(0.0000001)
        
    end
    
  
    % stop simulation if it blows up.
    if(norm(X)> 1e2)
      break
    end

end

%figure(3); hold on
%plot(timevec', (bottom_position(:,1).^2 + bottom_position(:,2).^2),'linewidth',2)
%plot(timevec', (bottom_position(:,1).^2 + bottom_position(:,2).^2),'linewidth',2)
%axis equal

% function which defines the force on the ground, existing on the plane X_3 = 0
function Fg = F_ground(X)
    % strength of the force exerted by the ground
    Sg = 800000000; 
    n=30;
    Fg = zeros(size(X));
    Fg(:,3) = (X(:,3) < 0).*(-Sg*X(:,3));
    %{
    if X(n+3,3)<0
        magnitude = Sg*(norm(X(n+3,:)));
    else
        magnitude = 0;
    end
    direction = (X(n+1,:)-X(n+3,:)) / norm(X(n+1,:)-X(n+3,:));
    Fg(n+3,:) = magnitude * direction;
    %}
    
    Fg(n+3,3) = 33*1000000;
   
end

% try to use a translation to fix the x and y coordinates of the
% bottom point
function X_prime = horizontal_translate(X)
    % the bottom point should have x and y coordinates 0 and 0
    % if its current position is (x, y, z), then we add a vector (-x, -y,
    % 0) to the whole system.
    n=30;
    trans_vec = X(n+3,:);
    trans_vec(3)=0;
    trans_vec = (-1).*trans_vec;
    X_prime = X + trans_vec;
end


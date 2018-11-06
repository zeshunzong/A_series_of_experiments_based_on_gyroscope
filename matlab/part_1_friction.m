
dim = 3;

% set time vector
dt = 1e-7;
end_time = 4e-1; 
timevec = 0:dt:end_time;

% set the number of nodes and links
n = 30;     % num of nodes on the wheel
num_nodes = n+3;
num_links = 2*n+2;

% jj&kk give the two nodes connected by a particular link
jj = zeros(num_links,1);
kk = zeros(num_links,1);

% set the position matrix & the velocity matrix
X = zeros(num_nodes, dim);
U = zeros(num_nodes, dim);

% mass of each node and in total
mass_vec = ones(num_nodes,1);
mass_total = 1*num_nodes;

angle = pi/12;      % stanting angle of the axis
half_len = 0.8;     % half length of the axis
rad_wheel = 1;      % radius of the wheel
friction_const = 5*10^10; % friction constant
grav_const = 1000000; % gravity constant
node_grav = zeros(num_nodes, dim);
node_grav(:,1) = 0;
node_grav(:,2) = 0;
node_grav(:,3) = -grav_const.*mass_vec;

% the rotation matrix is used to tilt the gyroscope into the starting
% position
rotation_along_y = [cos(angle), 0, sin(angle);0,1,0;-sin(angle), 0, cos(angle)];    

% set the starting position of each node
    % position of the nodes standing straight
    center = [0, 0, half_len];
    for k = 1:n
        theta = 2*pi*k/n;
        X(k,:) = center + rad_wheel*[cos(theta), sin(theta), 0];
    end
    X(n+1,:) = center;                    % center of mass
    X(n+2,:) = center + [0,0,half_len];   % top of the axis
    X(n+3,:) = center - [0,0,half_len];   % bottom of the axis

    % tilt the gyroscope
    for k = 1:n+3
        X(k,:) = (rotation_along_y * X(k,:)')';
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

Xcm = (sum((mass_vec.*X))./sum(mass_vec))';   % position of the center of mass
Ucm = (sum((mass_vec.*U))./sum(mass_vec))';   % velocity of the center of mass
Xtilda = X - Xcm';    % distance vector of each node to the center of mass

% initialize the angular momentum
L = (rotation_along_y * [0; 0; 100000]);

% this variable is used to compute the movement of the bottom node in each
% iteration and thus determine the friction
bottom_pos_ori = X(n+3, :);

for t = 1:length(timevec)

    % compute moment of inertia tensor
    I = zeros(dim, dim);
    for k = 1:num_nodes
        I = I + mass_vec(k).*((norm(Xtilda(k,:))^2).*eye(dim) - Xtilda(k,:)'*Xtilda(k,:) );
    end     

    % compute the angular velocity
    Omega = I\L;
    
    
    % normalize angular velocity if it is nonzero 
    % get the new distance vector Xtilda
    if(norm(Omega) > 100*eps)
         unit_Omega = Omega/norm(Omega);
         Omega_cross = [0 -Omega(3) Omega(2); Omega(3) 0 -Omega(1); -Omega(2) Omega(1) 0];
         P_Omega = unit_Omega*unit_Omega';
         Xtilda = (P_Omega*(Xtilda') + cos(norm(Omega)*dt).*(eye(dim) - P_Omega)*(Xtilda') + sin(norm(Omega)*dt).*(Omega_cross*(Xtilda'))./norm(Omega) )';
    end

    % compute net force and net torque   
        net_force = zeros(dim,1);
        net_torque = zeros(dim,1);
        
        % to simplift the algorithm, we cheat by assuming the normal force
        % equals the gravity and the net force is the normal force
        normal_force = zeros(num_nodes, dim);
        normal_force(n+3, :) = [0, 0, mass_total*grav_const];
        
        % compute the displacement of the bottom node and friction
        % we normalize it since we only need the direction
        bottom_mov = X(n+3, :)-bottom_pos_ori;
        if norm(bottom_mov)~=0
            bottom_mov = bottom_mov/norm(bottom_mov);
        end
        
        friction = -friction_const.*(bottom_mov);
        
        % construct a force vector for each node
        force_vec = node_grav+normal_force;
        force_vec(n+3, :) = force_vec(n+3, :)+friction;
        
        for l = 1:num_nodes
            net_force = net_force + force_vec(l, :)';
            net_torque = net_torque + cross(Xtilda(l,:)', force_vec(l, :)');
        end
        
        
        
    % update the position and velocity for the center of mass
    Ucm = Ucm + (dt/sum(mass_vec)).*net_force;
    Xcm = Xcm + dt.*Ucm;
    
    % update the angular momentum
    L = L + dt.*net_torque; 
    
    % record the position of the bottom node
    bottom_pos_ori = X(n+3, :);
    
    % update positions of individual masses
    X = Xtilda + Xcm';
    
    
    % plot
    mark = mod(t,1000);
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
    if(norm(X(:,:,1))> 1e2)
      break
    end

end


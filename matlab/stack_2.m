

% ----------------------------------SET UP----------------------------------
% set time vector
dt = 1e-4;
end_time = 2; 
timevec = 0:dt:end_time;

% set the number of nodes and links
num_gyros = 2;      % num of gyroscopes stacked together
n = 30;     % num of nodes on the wheel
num_nodes = n+3;
num_links = 2*n+2;

% jj&kk give the two nodes connected by a particular link
jj = zeros(num_links,1);
kk = zeros(num_links,1);

% set the position matrix & the velocity matrix
dim = 3;
X = zeros(num_nodes, dim, num_gyros);
U = zeros(num_nodes, dim, num_gyros);

% mass of each node and in total
mass_vec = ones(num_nodes, num_gyros);  % the i-th COLUMN is the mass vec of the i-th gyro
mass_total = sum(mass_vec, 1)';
% coefficients used to compute inter-gyro forces
    tension_const = 10^5;
    % Note that we will not multiply the dashpot constant by velocity of
    % the spring since to compute velocity we'll need to divide the
    % change of length by dt and that creats instability because dt is
    % small. Instead, we simply multiply it with the change of length
    % itself.
    dashpot_const = 100; 
    friction_const = 1; % friction constant
    normal_const = 10^5;

angle = pi/30;      % stanting angle of the axis
len = 1.8;          % length of the axis
half_len = 0.9;     % half length of the axis
rad_wheel = 1;      % radius of the wheel
grav_const = 10; % gravity constant
node_grav = zeros(num_nodes, dim, num_gyros);
for i = 1:num_gyros
    node_grav(:,1,i) = 0;
    node_grav(:,2,i) = 0;
    node_grav(:,3,i) = -grav_const.*mass_vec(:,i);
end

% the rotation matrix is used to tilt the gyroscope into the starting position
rotation_along_y = [cos(angle), 0, sin(angle);0,1,0;-sin(angle), 0, cos(angle)];    

% set the starting position of each node
    % the i-th row of the center matrix is the coordinate of the center of
    % the i-th gyroscope
    center = ones(num_gyros, dim);
    % position of the nodes standing straight
    for i=1:num_gyros
        center(i, :) = [0, 0, len*(i-1)+half_len];
        for k = 1:n
            theta = 2*pi*k/n;
            X(k, :, i) = center(i, :) + rad_wheel*[cos(theta), sin(theta), 0];
        end
        X(n+1, :, i) = center(i, :);                    % center of mass
        X(n+2, :, i) = center(i, :) + [0,0,half_len];   % top of the axis
        X(n+3, :, i) = center(i, :) - [0,0,half_len];   % bottom of the axis
    end
    % tilt the gyroscope 
    for i=1:num_gyros
        for k = 1:n+3
            X(k,:,i) = (rotation_along_y * X(k,:,i)')';
        end
    end
    % tilt all the gyroscopes back except the first one
    rotation_back = [cos(-angle), 0, sin(-angle);0,1,0;-sin(-angle), 0, cos(-angle)];
    for i=2:num_gyros
        for k = 1:n+3
            X(k,:,i) = (rotation_back * (X(k,:,i)-X(n+2,:,1))')'+X(n+2,:,1);
        end
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

% i-th ROW is the position of CM of the i-th gyro
Xcm = ones(num_gyros, dim);
% i-th ROW is the velocity of CM of the i-th gyro
Ucm = ones(num_gyros, dim);
% distance vector of each node to the center of mass (ROW vector of each page)
Xtilda = ones(num_nodes, dim, num_gyros);

for i = 1:num_gyros
    Xcm(i,:) = (sum(mass_vec(:,i).*X(:,:,i))./sum(mass_vec(:,i)));
    Ucm(i,:) = (sum(mass_vec(:,i).*U(:,:,i))./sum(mass_vec(:,i)));
    Xtilda(:,:,i) = X(:,:,i) - Xcm(i,:);    
end












% ----------------------------------SPIN----------------------------------
% some variables used in the iteration
L = ones(dim, num_gyros);   % the i-th COLUMN vector is the angular momentum of the i-th gyro
net_force = zeros(dim, num_gyros);   % the i-th COLUMN vector is the force on CM of the i-th gyro
net_torque = zeros(dim, num_gyros);   % the i-th COLUMN vector is the torque of the i-th gyro

% Spring's length and the change of length in unit time
spr_len = zeros(1, num_gyros-1);    % i-th entry is the spring connecting the i-th and (i+1)-th gyro
spr_len_change = zeros(1, num_gyros-1);

% initialize the angular momentum
L(:, 1) = rotation_along_y * [0; 0; 300];
L(:, 2) = [0; 0; -300];

% this variable is used to compute the movement of the bottom node in each
% iteration and thus determine the friction
bottom_pos_ori = X(n+3, :, 1);

%start video
v = VideoWriter('stack');
open(v);

% everything below is in this loop
for t = 1:length(timevec)
    %{
    First compute moment of inertia. Next compute w=(I^-1)L. The use w to 
    compute the new Xtilda.
    %}
        % compute moment of inertia tensor
        I = zeros (dim,dim,num_gyros);
        for i=1:num_gyros
            for k = 1:num_nodes
                I(:,:,i) = I(:,:,i) + mass_vec(k,i).*((norm(Xtilda(k,:,i))^2).*eye(dim) - Xtilda(k,:,i)'*Xtilda(k,:,i) );
            end     
        end
        % compute the angular velocity
        Omega = zeros(dim,num_gyros);
        for i=1:num_gyros
            Omega(:,i) = I(:,:,i)\L(:,i);
        end
        % normalize angular velocity if it is nonzero 
        % get the new distance vector Xtilda
        for i=1:num_gyros
            if(norm(Omega(:,i)) > 100*eps)
                 unit_Omega = Omega(:,i)/norm(Omega(:,i));
                 Omega_cross = [0 -Omega(3,i) Omega(2,i); Omega(3,i) 0 -Omega(1,i); -Omega(2,i) Omega(1,i) 0];
                 P_Omega = unit_Omega*unit_Omega';
                 Xtilda(:,:,i) = (P_Omega*(Xtilda(:,:,i)') + cos(norm(Omega(:,i))*dt).*(eye(dim) - P_Omega)*(Xtilda(:,:,i)') + sin(norm(Omega(:,i))*dt).*(Omega_cross*(Xtilda(:,:,i)'))./norm(Omega(:,i)) )';
            end
        end
    
    
    
    % compute net force and net torque
    for i=1:num_gyros   
        net_force = zeros(dim,1);
        net_torque = zeros(dim,1); 
        % construct a force vector for each node
        % we start with only gravity and gradually add other forces
        force_vec = node_grav(:,:,i);
            
        % Force applied to the top node
        % if it's the top gyro, it's 0
        % otherwise, it's generated by a spring
        force_top = zeros(1, dim);
        if i ~= num_gyros
            top_str_len = X(n+3,:,i+1)-X(n+2,:,i);
            force_top = (tension_const-dashpot_const*spr_len_change(i)) * top_str_len;
        end
        force_vec(n+2,:) = force_vec(n+2,:)+force_top;
            
        % Force applied to the bottom node
        % if it's the bottom gyro, it's the normal force
        % otherwise, it's generated by a spring
        force_bot = zeros(1, dim);
        if i==1
            % we first compute the normal force
            normal_force = zeros(1, dim);
            if X(n+3,3,i)<0
                normal_force = -normal_const*[0,0,X(n+3,3,1)];
            end
            force_bot = force_bot+normal_force;   
            % next, we use the displacement of the bottom node
            % to compute friction
            bottom_mov = X(n+3,:,1)-bottom_pos_ori;
            if norm(bottom_mov)~=0
                bottom_mov = bottom_mov/norm(bottom_mov);
            end    
            %{ 
            friction is porportional to the normal force 
            note that we project the displacement of the bottom node onto
            the XY-plane in order to make friction horizontal
            %} 
            friction_magnitude = -friction_const*norm(normal_force);
            friction_direction = dot(bottom_mov, [1,0,0]).*[1,0,0] ...
                +dot(bottom_mov, [0,1,0]).*[0,1,0];
            friction = friction_magnitude.*friction_direction;
            force_bot = force_bot+friction;
            % update the new position of the bottom node
            bottom_pos_ori = X(n+3,:,1);
        else
            bot_str_len = X(n+2,:,i-1)-X(n+3,:,i);
            force_bot = (tension_const-dashpot_const*spr_len_change(i-1)) * bot_str_len;
        end
        force_vec(n+3,:) = force_vec(n+3,:)+force_bot;

        % compute the total force and torque acting on each gyroscope
        for l = 1:num_nodes
            net_force = net_force + force_vec(l, :)';
            net_torque = net_torque + cross(Xtilda(l,:,i)', force_vec(l, :)');
        end
        % update the position and velocity for the center of mass
        Ucm(i,:) = Ucm(i,:) + (dt/mass_total(i)).*net_force';
        Xcm(i,:) = Xcm(i,:) + dt.*Ucm(i,:);
        % update the angular momentum
        L(:,i) = L(:,i) + dt.*net_torque; 
    end
    
    
    
    % update positions of individual gyros
    for i=1:num_gyros
        X(:,:,i) = Xtilda(:,:,i) + Xcm(i,:);
    end
    % update the spring lengths
    spr_len_new = zeros(1, num_gyros-1);    
    for j=1:(num_gyros-1)
        spr_len_new(j) = norm(X(n+2,:,j)-X(n+3,:,j+1));
    end
    spr_len_change = spr_len_new - spr_len;
    spr_len = spr_len_new;
    
    
    
    % graph
    mark = mod(t,100);
    for i=1:num_gyros
        if mark < 1
            figure(1);
            x = zeros(num_links*num_gyros, 2);
            y = zeros(num_links*num_gyros, 2);
            z = zeros(num_links*num_gyros, 2);
            for j=1:num_gyros
                x(((j-1)*num_links+1):(j*num_links),:) = [X(jj,1,j) X(kk,1,j)];
                y(((j-1)*num_links+1):(j*num_links),:) = [X(jj,2,j) X(kk,2,j)];
                z(((j-1)*num_links+1):(j*num_links),:) = [X(jj,3,j) X(kk,3,j)];
            end
            plot3(x',y',z','linewidth',3)

            axis equal
            view([0,5])
            xlim([-3 3])
            ylim([-1 1])
            zlim([-2 6])
            pause(0.00000001)

            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end
end

close(v)

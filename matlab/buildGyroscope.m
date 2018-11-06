function [X, jj, kk, U, M, Xcm, Ucm, Xtwiddle, L] = buildGyroscope(n, unitvec1, Rbase, gyropos)
    % let us attach a gyroscope to this track, we do so by first build a
    % gyroscope centered at the origin, then change its orientation and shift 
    % its center of mass
    % this gyroscope should have a radius of Rbase


    % number of nodes on the circumference of the wheel
    % n = 30;

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
    X(n+3,:) = center - [0,0,Rbase]; % it should be that the bottom is at the origin

    % STEP TWO: rotate the structure
    % we want to rotate the gyroscope so that its rotating axis is parallel to
    % the direction of unitvec1.
    % we do so by (change of coordinates) -> (rotation) -> (change back)

    avec = [unitvec1(1), unitvec1(2), 0]; % projection of unitvec1 on xOy plane
    avec = avec ./ norm(avec); % normalize it, this is the new x direction
    cvec = [0 0 1]; % we still use the original z direction 
    bvec = cross(cvec, avec); % construct the other direction by cross prod
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
    
    % assume each node has mass 1
    M = ones(num_nodes,1); 
    % initialize the position of the center of mass
    Xcm = (sum((M.*X))./sum(M))';
    % initialize velocity center of mass
    Ucm = (sum((M.*U))./sum(M))';

    % initialize Xtwiddle
    %Xtwiddle = zeros(num_nodes,NDIM);
    Xtwiddle = X - Xcm';
    
    % also rotate the angular momentom so that it is parallel with the
    % rotating axis
    L = action * [0; 0; 100000];
end
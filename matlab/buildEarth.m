function [Xearth, jje, kke] = buildEarth(Rearth, nseg, ncir)
    % we first need to build an earth
    % we do so by first using spherical coordinates (R, theta, phi)
    %Rearth = 100;
    % we divide each latitude line into nseg segments
    %nseg = 12;
    thetas = 0 : 2*pi/nseg : 2*pi - 0.001;
    % we cut the ball horiontally to get ncir circles
    %ncir = 7;
    phis = pi/(ncir+1) : pi/(ncir+1) : pi-0.001;
    % so in total 5 circles, each circle has 6 points on it, 30 nodes
    % the nodes are indexed from the bottom to the top, one circle after the
    % other

    Xearth = zeros(ncir*nseg, 3);
    for i = 1:ncir
        % loop through the circles
        for j = 1:nseg
            % loop through the points on the circles
            % x = r*cos(theta)*sin(phi); y = r*sin(theta)*sin(phi); z =
            % r*cos(phi).
            Xearth(nseg*(i-1)+j,:) = [Rearth*cos(thetas(j))*sin(phis(i)), Rearth*sin(thetas(j))*sin(phis(i)), Rearth*cos(phis(i))];
        end
    end
    % how many links are there? on each circle there are nseg links, so nseg*ncir links here.


    jje = zeros(ncir*nseg,1);
    kke = zeros(ncir*nseg,1);
    % draw latitudes
    % this should fill in the first nseg*ncir links
    for i = 1:ncir 
        jje((i-1)*nseg+1 : i*nseg) = (i-1)*nseg+1 : i*nseg;
        kke((i-1)*nseg+1 : i*nseg) = [(i-1)*nseg+2 : i*nseg, (i-1)*nseg+1];
    end
end
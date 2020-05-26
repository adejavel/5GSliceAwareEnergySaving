tic;

% COPY HERE THE CONTENT OF ONE OF THE PARAMETERS FILE DEPENDING ON THE NEEDED SLICE
% eMBB: params_emmb.txt
% URLLC: params_urllc.txt
% MIoT: params_miot.txt
% All slices: params_all.txt
% Or create your own file by copying one of the previous file


% Computing slices BWP
capacity_per_surface = [];
total_capacity_per_surface = 0;
for s=1:size(capacity_per_users,2)
    % we get capacity per user requirement
    slice_capacity = capacity_per_users(s);
    % We get slice UEs density
    slice_density = densities(s);
    % We get slice UEs connected mode time proportion
    slice_active_proportion = active_proportions(s);
    slice_per_surface_capacity = slice_density * slice_active_proportion * slice_capacity;
    
    capacity_per_surface = [capacity_per_surface slice_per_surface_capacity];
    total_capacity_per_surface = total_capacity_per_surface + slice_per_surface_capacity;
end
bwp_limit = [];
numerologies = [];
number_rbs = [];
for s=1:size(capacity_per_users,2)
    % we get capacity per user requirement
    slice_capacity = capacity_per_users(s);
    % We get slice UEs density
    slice_density = densities(s);
    % We get slice UEs connected mode time proportion
    slice_active_proportion = active_proportions(s);
    % We get slice latency criticity
    slice_latency_criticity = latency_criticity(s);
    % Get per surface capacity
    slice_per_surface_capacity = capacity_per_surface(s);
    
    % compute the proportion of bandwidth that can be allocated to the
    % slice
    allowed_bandwidth = bandwidth * slice_per_surface_capacity / total_capacity_per_surface;
    bwp_limit = [bwp_limit allowed_bandwidth];
    
    % determining slice BWP SCS
    if slice_latency_criticity==0
        slice_numerology = 15e3;
    elseif slice_latency_criticity ==1
        slice_numerology = 30e3;
    else
        slice_numerology = 60e3;
    end
    numerologies = [numerologies slice_numerology];
    rb_size = 12 * slice_numerology;
    % determining slice BWP number of RBs
    slice_number_rb = fix(allowed_bandwidth/rb_size);
    if slice_number_rb < 1
        slice_number_rb = 1
    end
    number_rbs = [number_rbs slice_number_rb];
end


slices_params = [densities ; active_proportions ;capacity_per_users;numerologies;number_rbs ];

max_connectivity = max(k_connectivity);
reference_distance = 0.05*sqrt(1/lambda);

rng('default')
num_sim_done = 0;

while num_sim_done<num_sim
    parfor sim = 1:num_sim-num_sim_done
        rng(round(second(now))+sim);
        % nbAntennas is the number of antennas which follows a Poisson Point
        % Process

        nbAntennas = random('Poisson',mean_antennas)
        if nbAntennas < 1
            nbAntennas = 1
        end



        % Variables used to plot r_max and radius depending on the power
        radii_1 = [];
        powers_1 = [];

        % Border antennas and border radii are containing the antennas that are
        % enclosing the studied zone. Those antennas are taken in account for the
        % Homology part but not for the SINR computation


        border_antennas = [];
        border_radii = [];

        for i=0:100:border
            border_antennas = [border_antennas; [i 0]];
            border_radii = [border_radii; 60];
            border_antennas = [border_antennas; [i border]];
            border_radii = [border_radii; 60];
            if i ~= 0 && i ~= border
                border_antennas = [border_antennas; [0 i]];
                border_radii = [border_radii; 60];
                border_antennas = [border_antennas; [border i]];
                border_radii = [border_radii; 60];
            end
        end



        % Generating the antennas positions using the number of antennas
        % 15 offset avoid having lot of antennas in the border

        antennas = randi([30 border-30],nbAntennas,2);
        valid_distrib = true;
        for a1=1:size(antennas,1)
            for a2=a1+1:size(antennas,1)
                distance = sqrt((antennas(a1,1)-antennas(a2,1))^2+(antennas(a1,2)-antennas(a2,2))^2)
                if distance < reference_distance
                    valid_distrib = false;
                    break;
                end
            end
        end
        


        if valid_distrib
            num_sim_done = num_sim_done +1; 
            % Emission power contains the emission power for all the antennas (i-th element of this arradii represents the emission power of the i-th antenna)
            % The antennas are currently initialised at 10W
            emission_power = zeros(1,size(antennas,1)) + base_emission_power;
            % First radii computation to initialize our algorithm
            radii = compute_radii(antennas,emission_power,radio_params,slices_params);

            % total_antennas and total_radii contain the studied antennas and the border
            % antennas 
            total_antennas = [ antennas; border_antennas ];
            total_radii = [radii; border_radii];

            % First holes (Betti number) computation 
            % From now on, we want our algorithm to maintain this number of holes
            % Holes is an arradius containing the betie numbers (=> k-dimension coverage hole)
            [holes,k_simplex] = compute_holes(total_antennas,total_radii,max_connectivity);

            % Pruning all the alone antennas
            simplexes1 = k_simplex{1,2};
            new_antennas = [];
            seen_simplexes = [];
            for r=1:size(simplexes1,1)
                simplex01 = simplexes1(r,1)
                simplex02 = simplexes1(r,2)
                if simplex01 <= nbAntennas && ~ ismember(simplex01,seen_simplexes)
                    new_antennas = [new_antennas;antennas(simplex01,:)];
                    seen_simplexes = [seen_simplexes simplex01];
                end
                if simplex02 <= nbAntennas && ~ ismember(simplex02,seen_simplexes)
                    new_antennas = [new_antennas;antennas(simplex02,:)];
                    seen_simplexes = [seen_simplexes simplex02];
                end
            end

            antennas = new_antennas;

            emission_power = zeros(1,size(antennas,1)) + base_emission_power;
            radii = compute_radii(antennas,emission_power,radio_params,slices_params);

            total_antennas = [ antennas; border_antennas ];
            total_radii = [radii; border_radii];

            [holes,k_simplex] = compute_holes(total_antennas,total_radii,max_connectivity);




            beginning_radii = total_radii;

            new_holes = holes;

            % Beginning simulated annealing algorithm
            % t0 is the inital temperature
            % alpha ^ k represents the cooling factor

            t0 = (-DP) / log(p_init);
            alpha = nthroot( (-DP) / ( t0 * log(p_end) ) , K );



            st_count = 0;

            % for each step
            for k=0:K
                % compute step temperature tk
                tk = power( alpha , k ) * t0
                % deduct the probability for a power increase to be accepted
                prob_threshold = exp( - DP / tk );
                % for each sub-step
                for l = 1 : L
                    st_count = st_count + 1;
                    disp("###############")
                    disp("SIMU " +sim+"/"+num_sim+" // STEP "+ st_count + " / "+ (K+1)*L )
                    disp("###############")
                    % choose an antenna
                    cell = randi([1 size(antennas,1)]);
                    % choose 0 or 1 (0 => power decrease, 1 => power increase) 
                    sign = randi([0 1]);
                    % Power decrease
                    if sign == 0
                        % first we decrease the power of the selected cell by dp 
                        emission_power(cell) = emission_power(cell) - dp;
                        % We compute new radii for all our antennas 
                        radii = compute_radii(antennas,emission_power,radio_params,slices_params);
                        % For holes computation, we need our borders so we add border
                        % radii
                        total_radii = [radii; border_radii];
                        % we compute our new holes
                        [init_new_holes,k_simplex] = compute_holes(total_antennas,total_radii,max_connectivity)
                        % adjust_holes is a function that adjust arradius sizes ([1 0] and [2 3 4] would be returned as [1 0 0] and [2 3 4])
                        [holes,new_holes] = adjust_holes(holes,init_new_holes);
                        % Check k-connectivity

                        max_index = min(size(init_new_holes,2),max_connectivity+1);

                        if ~ all(new_holes(1:max_index) <= holes(1:max_index)) || ~ all(emission_power >= min_power)
                            emission_power(cell) = emission_power(cell) + dp;
                        end
                    % power increase 
                    else
                        % we take a random number between 0 and 1.
                        % the chance for this random number to be smaller than
                        % prob_threshold represents the increase validation probability
                        % we also check that our emission power will not be higher than
                        % 40W
                        x = rand;
                        if x < prob_threshold && emission_power(cell) < 40 - DP
                            % if our increase is accepted, we increase the cell power
                            % by DP
                            emission_power(cell) = emission_power(cell) + DP;
                        end
                    end
                    radii_1 = [radii_1 radii(1)];
                    powers_1 = [powers_1 emission_power(1)];
                end
            end



            % final steps
            % we compute our final radii and we plots different results
            radii = compute_radii(antennas,emission_power,radio_params,slices_params);
            total_radii = [radii; border_radii];

            rjson = [];
            rjson.emission_power = emission_power;
            rjson.beginning_radii = beginning_radii;
            rjson.end_radii = total_radii;
            rjson.holes = holes;
            rjson.k_simplex = k_simplex;
            rjson.nbAntennas = nbAntennas;
            rjson.total_antennas = total_antennas;


            global_result = [global_result rjson];
        end

    end
end



pjson.k_connectivity = k_connectivity;
pjson.K = K;
pjson.L = L;
pjson.dp = dp;
pjson.DP = DP;
pjson.p_init = p_init;
pjson.p_end = p_end;
pjson.lte_threshold = lte_threshold;
pjson.band = band;
pjson.gamma = gamma;
pjson.c = c;
pjson.bandwidth = bandwidth;
pjson.nb_circles = nb_circles;
pjson.densities = densities;
pjson.active_proportions = active_proportions;
pjson.capacity_per_users = capacity_per_users;
pjson.latency_criticity = latency_criticity;
pjson.numerologies = numerologies;
pjson.number_rbs = number_rbs;
pjson.lambda = lambda;
pjson.base_emission_power = base_emission_power;


toc;



% This function is fitting holes and new_holes sizes
function [holes, new_holes] = adjust_holes(holes,new_holes)
    while size(holes,1) < size(new_holes,1)
        holes = [holes ; 0];
    end
    while size(holes,1) > size(new_holes,1)
        new_holes = [new_holes; 0];
    end
    
end


% This function computes cell capacity radius for each antenna, given the
% emission power
function radii = compute_radii(antennas, emission_power,radio_params,slices_params)

    % antenna_x and antenna_y contain all the X and Y coordinates of the
    % antennas
    antennas_x = antennas(:,1);
    antennas_y = antennas(:,2);

    % radius arradius will contain all the radii of the cells
    radii = [];
    
    lte_threshold = radio_params(1);
    band = radio_params(2);
    gamma = radio_params(3);
    c = radio_params(4);
    bandwidth = radio_params(5);

    densities = slices_params(1,:);
    active_proportions = slices_params(2,:);
    capacity_per_users = slices_params(3,:);
    numerologies = slices_params(4,:);
    number_rbs = slices_params(5,:);

    nb_circles = radio_params(6);
    
    % We first determine our reference distance depending on the used band 
    % The basic equation is
    % r0 = 2 * band * D ^ 2 / c; 
    % with D the antenna size but we consider an antenna size of 1m (=> D=1)
    r0 = 1; 
    
    disp("beginning radius calcul");
    % for each cell
    parfor i=1:size(antennas,1)
        % if the antenna is ON
        if emission_power(i)>0
            % We first get the antenna coordinates x_cell and y_cell
            x_cell = antennas_x(i);
            y_cell = antennas_y(i);
           
            % We initalize our current_radius variable to the minimum
            % refenrece distance
            current_radius =  r0;
            
            % r_max is the maximum radius for which the UE will receive the
            % signal with a power of lte_threshold. Beyond this r_max,
            % reception power will not be big enough for signal decoding
            % To determine this radius, we use the path loss model
            r_max = 2 * band / ( c * nthroot((( lte_threshold / emission_power(i) ) * ( 8 * pi ) ^ 2 * ( band / c ) ^ 4), gamma));
            
            % We are now building artificial rings (we build nb_circles rings) inside the circle
            % defined by r_max. Each ring will be studied independently 
            
            % As the capacity decreases exponentially over distance, we want our
            % rings radius to increase exponentially so the capacity is linear
            % over rings
            % We build an exponential rule called exp such as 
            % exp ^ (nb_circles) = r_max
            exp = nthroot(r_max,nb_circles);
            % For each ring
            for z=1:nb_circles
                % k is the ring radius
                k = exp^z;
                % We only consider radii higher than the reference distance
                if k > r0
                    % We first determine the reception power at the ring
                    % edge using the path loss model
                    pr = emission_power(i) * power( 2 * band / ( k * c ), gamma) * power( c ^ 2 / ( 8 * pi * band ^ 2 ) , 2);
                    % Then we compute interferences (=> sum of reception powers for all other cells)
                    % We initialize interferences to 0
                    
                    interferences = 0;
                    % for all neighbor cells
                    for j=1:size(antennas,1)
                        % (except the studied cell)
                        if j ~= i
                            % get neighbor cell coordinates
                            x_neigh = antennas_x(j);
                            y_neigh = antennas_y(j);
                            % compute distance between the studied cell and
                            % the neighbor cell
                            Dk = sqrt(power(x_cell - x_neigh, 2) + power(y_cell - y_neigh, 2));
                            % dk is the distance for which neighbor cell interference
                            % will be the higher (i.e. on the edge of the ring)
                            dk = abs(Dk - k);
                            % then we compute the reception power from
                            % this neighbor cell using path loss modle and
                            % dk distance
                            new_int = emission_power(j) * power( 2 * band / ( dk * c ), gamma) * power( c ^ 2 / ( 8 * pi * band ^ 2 ) , 2);
                            % and finally we sum it to other interferences
                            interferences = interferences + new_int;
                        end                
                    end
                   
                    is_ring_validated = true;
                    % For each slice
                    for s=1:size(capacity_per_users,2)
                        
                        % we get capacity per user requirement
                        slice_capacity = capacity_per_users(s);
                        % We get slice UEs density
                        slice_density = densities(s);
                        % We get slice UEs connected mode time proportion
                        slice_active_proportion = active_proportions(s);
                        slice_numerology = numerologies(s);
                        slice_number_rb = number_rbs(s);
                        % and finally we add to capacity_required the
                        % capacity required by the studied slice
                        capacity_required = slice_capacity * slice_active_proportion * slice_density * pi * k ^ 2;
                        % offered capacity
                        slice_bandwidth = slice_numerology * 12 * slice_number_rb;
                        noise = 1.38066e-23 * 290 * slice_bandwidth;
                        % We are now able to compute the SINR 
                        sinr = pr / (interferences+noise);
                        offered_capacity = slice_bandwidth * log2(1 + sinr);
                        if offered_capacity < capacity_required
                            is_ring_validated = false;
                            break;
                        end
                        
                    end
                    % Then, if our total_capacity offered by our ring is
                    % higher than the capacity_required by slices, our ring
                    % will be validated by saving its radius and we continue
                    if is_ring_validated
                        current_radius = k;  
                    % If ring is not validated, we stop our loop
                    else
                        break;
                    end
                end
            end
            % We finally add to cell radii the current_radius variable which
            % represents the radius of the bigger ring that serves slices
            % requirements
            radii = [radii; current_radius];
        % if the antenna is OFF
        else
            radii = [radii; 0];
        end
    end
end


% compute_holes function computes Betti numbers for our topology
function [holes,k_simplex] = compute_holes(antennas, radii,max_connectivity)
    % antenna_x and antenna_y contain all the X and Y coordinates of the
    % antennas
    antennas_x = antennas(:,1);
    antennas_y = antennas(:,2);

    % We first determine k-simplices
    
    % Zero-simplices are all the antennas
    zero_simplex = [];
    antennas_size = size(antennas);
    zero_simplex = 1:antennas_size(1);
    zero_simplex = zero_simplex.';
    
    % k_simplex is the object containing our k-simplex (k_simplex{i} contains (i-1)-simplices)
    k_simplex={};
    k_simplex{1}=zero_simplex;

    % 1-simplices computation
    % 1-simplices are pair of 0-simplices (antennas) whose coverage zone
    % intersection is not null
    % We first get all pairs of antennas
    combinations = all_combinations(zero_simplex,2);
    combinations_elements = combinations(:,:);
    % Then we get radii for each pair of elements
    combi_radii = radii(combinations_elements);
    % we compute coordinates for our pairs
    combinations_coord = [antennas_x(combinations_elements) antennas_y(combinations_elements)];
    % Then we create a condition for pair selection
    % The distance between the two antennas of the pair must be smaller
    % that the sum of the coverage radius of each antenna
    connections = sqrt(power(combinations_coord(:,1)-combinations_coord(:,2),2)+power(combinations_coord(:,3)-combinations_coord(:,4),2)) < (combi_radii(:,1)+combi_radii(:,2));
    % We get all pairs that fit this condition
    one_simplex = combinations_elements(connections,:);
    % And finally we store our 1-simplices
    k_simplex{2} = one_simplex;


    % Then we compute all other simplices with the same algorithm
    % To build k-simplices we merge (k+1) (k-1)-simplices
    k=2;
    % while we are able to build (k+1)-simplices
    while size(k_simplex,2)>=k && size(k_simplex{k},1)>=k+1 && k <= max_connectivity + 1
        good_size_arrays = [];
        count_array=[];
        first = true;
        for p=1:size(k_simplex{k},1)
            for n=1:size(k_simplex{k}(p+1:size(k_simplex{k},1),:),1)
                n_array = unique([k_simplex{k}(p,:) k_simplex{k}(n+p,:)]);
                if size(n_array,2)==k+1
                    %good_size_arrays = [good_size_arrays; n_array];
                    
                    if first
                        good_size_arrays = [good_size_arrays; n_array];
                        count_array = [count_array 1];
                        first = false;
                    else
                        [res, posi] = ismember(n_array,good_size_arrays,'rows');
                        if res
                            count_array(posi) = count_array(posi)+1;
                        else
                            good_size_arrays = [good_size_arrays; n_array];
                            count_array = [count_array 1];
                        end
                    end
                    
                end
            end

        end
        comb_numb = factorial(k+1)/(factorial(2)*factorial(k-1));
        k_simplex{k+1} = good_size_arrays(count_array==comb_numb,:);
        k = k+1;
    end
    % Then we compute matrices that represent bases for edge operator
    matrices = {};
    parfor j=1:size(k_simplex,2)-1
        matrices{j} = zeros([size(k_simplex{j},1) size(k_simplex{j+1},1)]);
        for i=1:size(k_simplex{j+1},1)
            simplex_element = k_simplex{j+1}(i,:);
            for k=1:j+1
                temp_simplex = simplex_element;
                temp_simplex(k) = [];
                [result, pos] = ismember(temp_simplex,k_simplex{j},'rows');
                matrices{j}(pos,i) = power(-1,k+1);
            end
        end
    end

    % And finally we compute Betti numbers stored in beta variable
    beta=[];
    beta = [beta; (size(k_simplex{1},1) - rank(matrices{1}))];
    for i=2:size(matrices,2)
        beta = [beta; (size(k_simplex{i},1) - rank(matrices{i-1}) - rank(matrices{i}))];
    end
    holes = beta;
end

% This functions computes all combinations of source
function array = all_combinations(source,combination_size)
    combinations = combntns(1:size(source,1),combination_size);
    array = zeros([size(combinations,1) combination_size size(source,2)]);
    parfor i = 1:size(combinations,1)
        for k = 1:combination_size
            array(i,k,:)=source(combinations(i,k),:);
        end     
    end
end



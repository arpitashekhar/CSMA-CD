function projectdemo3()
    clear all;
    tic
    
    % NODE Id's
     A = 1; B = 2; C = 3; D = 4;
    nodes = [A, B, C, D];
    % Elist is an (n by 7) matrix where n is the number of packets
    % waiting in the queue of respective nodes. We define a packet to be a
    % row in the elist and is of the form
    % [SRC DEST GENTIME TXTIME RXTIME CURTIME COLLISIONS]
    elist1 = [];    % for bus 1 connecting Node A and B
    elist2 = [];    % for bus 2 connecting Node C and D
    
    % Columns in the event list
    SRC = 1;      % Id of the source node
    DEST = 2;     % Id of desstination node
    GENTIME = 3;  % Time at which the packet is created at the source node
    TXTIME = 4;   % Time at which the packet is transmitted
    RXTIME = 5;   % Time at which the packet is received at the dest node
    % The current simulation time w.r.t the packet. You can think of
    % CURTIME to be the time at which the simulator had last visited this
    % packet.
    CURTIME = 6;  
    COLLISIONS = 7; % Number of collisions incurred by 
                    % this packet due to channel contention
    
    % Reference:
    % http://en.wikipedia.org/wiki/Discrete_event_simulation#Clock 
    % Clock:
    % The simulation must keep track of the current simulation time,
    % whatever measurement units are suitable for the system being modeled.
    % In discrete-event simulations, as opposed to real time simulations,
    % time ‘hops’ because events are instantaneous – the clock skips to the
    % next event start time as the simulation proceeds. 
    CLOCK = 0;
 
    TOTALSIM = 30*10^3; % Total simulation time
    lambda = .5;
    frameslot = 500; % frame slot time (usec) <-- changed from 500us
    td = 80;% transmission delay on BUS (usec)
    pd = 10; % propagation delay on BUS
    tdelay = td + pd; % total delay incurred during a pkt transmission
	tbackoff = frameslot; % time slot (usec) for backoff algorithm
    maxbackoff = 3; % maximum backoff time is 2^3 frame slot
	
    % the time at which the last packet was generated at node A, B, C and D resp.
    GENTIMECURSOR = [0 0 0 0];
        
    % create packets for all the nodes
    [x, y] = size(nodes);
     for n = 1:y
         createpacket(n);
     end
    
    if(size(elist1, 1) == 0 && size(elist2, 1) == 0)
        disp('No packets to simulate');
        return;
    end
    
    % Collect the statistics in this array.
    SIMRESULT = [];
        
    while(1)
        
        % Update the clock.
        updateclock(); 
        
        % Find the source node of the packet randomly from the first row of elist1 and elist2.
        l1 = elist1(1,SRC); l2 = elist2(1,SRC);
        l = horzcat(l1, l2);
        randomindex = randi(length(l));
        src = l(randomindex);
		
        % Randomly generating destination nodes and adding to the elist.
        destnodes = nodes;
        destnodes(src) = [];
        randomindex = randi(length(destnodes));
        dst = destnodes(randomindex);
        
        % if source node is connected to bus 1, update elist1; otherwise
        % update elist2
        if src == 1 || src == 2
            bus1 = true; bus2 = false;      % bus 1 will be used for transmission
            elist1(1,DEST) = dst;
        else
            bus1 = false; bus2 = true;      % bus 2 will be used for transmission
            elist2(1,DEST) = dst;
        end
        
        % route the packets if the source and desination are connected to
        % different destinations
        if ((src == 1 || src == 2) && (dst == 3 || dst == 4)) || ((src == 3 || src == 4) && (dst == 1 || dst == 2))
            routing = true;
            routingpackets(src, dst);
        else
            routing = false;
        end
        
        timediff1 = elist1(2,CURTIME) - elist1(1,CURTIME);
        timediff2 = elist2(2,CURTIME) - elist2(1,CURTIME);
        
		if(timediff1 > pd && timediff2 > pd)
			%% No collision case
            % Set the tx time. This should be the time when the packet
            % is transmitted for the first time.
            
            if ~routing
                tdelay = 90;        % resetting tdelay to 90 when there is no routing, i.e. only bus is used for transmission
            end
            
            % if the packet is transmitted from elist1, update the
            % transmission and received time for the packet in elist1
            if bus1
                if elist1(1,TXTIME) == 0
                    elist1(1,TXTIME) = elist1(1,CURTIME);
                end
                % Set the rx time.
                elist1(1,RXTIME) = elist1(1,CURTIME) + tdelay;
            end
            
            % if the packet is transmitted from elist2, update the
            % transmission and received time for the packet in elist2
            if bus2
                if elist2(1,TXTIME) == 0
                    elist2(1,TXTIME) = elist2(1,CURTIME);
                end
                % Set the rx time.
                elist2(1,RXTIME) = elist2(1,CURTIME) + tdelay;
            end

            updatesimlist();
			
            createpacket(src);
			
			% add tdelay to CLOCK. Check delaypkts() for more details.
            delaypkts(tdelay);
       else
            backoffoncollision();
       end
		
        if min(GENTIMECURSOR) > TOTALSIM
            disp('Completed!');
            calcstat();
            break;
        end
    end
      % calculating statistics for all the packets that are simulated 
	   function calcstat()
        % nodename = ['A', 'B', 'C', 'D']; 
        figure(1);
        figure(2);
        plotcount = 0;
        f = [];
        for i = 1:4		% source
			source = nodes(i);
            var3 = 0;
            queuedelay = 0;
            accessdelay = 0;
			for j = 1:4		% destination
				if j ~= i
                    plotcount = plotcount + 1;
                    % calculating packet statistics
					destination = nodes(j);                    
					var1 = SIMRESULT(SIMRESULT(:,SRC) == source,:);
                    var1 = var1(var1(:,DEST) == destination,:);
                    
                    % number of packets sent from node i to node j
                    var2 = length(var1);
                    
                    % queuing delay from node i to node j
                    var3 = var1(:, RXTIME) - var1(:, GENTIME);
                    figure(1);
                    plottitle = strcat('Queue delay from node ', num2str(i), ' to ', num2str(j));
                    queuedelay = sum(var3);     % Sum of queuing delays at every node
                    subplot(4,3, plotcount);
                    plot(1:length(var3),var3);
                    xlabel('Packet sequence #');
                    ylabel('Delay in \mu sec');
                    title(plottitle);
                    
                    % Access delay
                    var4 = var1(:, RXTIME) - var1(:, TXTIME);
                    figure(2);
                    accessdelay = sum(var4);
                    plottitle = strcat('Access delay from node ', num2str(i), ' to ', num2str(j));
                    subplot(4,3, plotcount);
                    plot(1:length(var4),var4);
                    xlabel('Packet sequence #');
                    ylabel('Delay in \mu sec');
                    title(plottitle);
                    
                    var5 = var1(:, COLLISIONS);
                    noofcollisions = sum(var5);
                    % end to end throughput for all pair of nodes
                    meanendtoend = mean(var3);
                    
                    % Average end to end throughput
                    % Avgtha=((1000*8)/meanendtoenda)*10^6;  % bits/sec
                    avgthroughput = ((1000*8)/meanendtoend)*10^6;  % bits/sec
                    f = [f; i j avgthroughput var2 noofcollisions];
               end
            end
            % total queuing delay at node A, B, C and D
            queuedelay = queuedelay - 3*tdelay;
        end
        disp(array2table(f, 'VariableNames',{'Source','Destination','Throughput', 'NoOfPacketsSent', 'NoOfCollisions'}));
    end
	
    % The clock slips to the RXTIME i.e., add delay time to CLOCK.
    function delaypkts(delay)
        CLOCK = CLOCK + delay;
        % It might so happen that the new packet at the SRC node might have
        % been generated when the previous packet was in flight. This new
        % packet cannot be transmitted immediately and hence has to wait
        % till the previos packet has reached the destination. 
        
        % list will have the row number of elist whose CURTIME field value
        % is less than CLOCK. Remember CLOCK is now the RXTIME. 
        list1 = find(((elist1(:,CURTIME)-CLOCK) < 0));
        list2 = find(((elist2(:,CURTIME)-CLOCK) < 0));
        % Set the CURTIME field of all the rows in list to CLOCK.
        elist1(list1,CURTIME) = CLOCK;
        elist2(list2,CURTIME) = CLOCK;
    end
   
    function updateclock()
        % SORTROWS(elist,CURTIME) sorts the rows of elist in ascending
        % order for the column CURTIME.
        elist1 = sortrows(elist1, CURTIME);
        elist2 = sortrows(elist2, CURTIME);
        
        % Set the clock to the CURTIME of the packet in the first row of
        % the elist since this is the packet that contends first for the
        % channel. 
        CLOCK = min(elist1(1,CURTIME), elist2(1,CURTIME));
    end
    
    % If the number of arrivals in any given time interval [0,t] follows
    % the Poisson distribution, with mean = \lambda \cdot t, then the lengths of the
    % inter-arrival times follow the Exponential distribution, with mean
    % 1/\lambda.
    function pkt = createpacket(nodeid)
        % Find the inter-arrival time.
        interarvtime = round(frameslot*exprnd(1/lambda,1,1));
        
        % Find the birth time.
        GENTIMECURSOR(nodeid) = GENTIMECURSOR(nodeid) + interarvtime;
        
        % Create the packet. Unknown fields are set to 0.
        % [SRC = nodeid DEST = 0 GENTIME = birthtime TXTIME = 0, RXTIME = 0
        % CURTIME = birthtime COLLISIONS = 0]
        pkt = [nodeid 0 GENTIMECURSOR(nodeid) 0 0 GENTIMECURSOR(nodeid) 0];
        
        % Enqueue to the event list.
        if(nodeid == 1 || nodeid == 2)
            elist1 = [elist1; pkt];
        else
            elist2 = [elist2; pkt];
        end
    end
	
	% Return the CURTIME of node
	function t = getcurtime(node, elist)
        % idx will have the row number of elist whose CURTIME field value
        % is less than DELAYTIME.
		idx = find(elist(:,SRC) == node, 1, 'first');
        t = elist(idx, CURTIME);
    end
	
    function delaynodepkts(node, delay)
        if bus1
            DELAYTIME = getcurtime(node, elist1) + delay;
            list = find(elist1(:,CURTIME)-DELAYTIME < 0 & elist1(:,SRC)==node); 
            % Set the CURTIME field of all the rows in list to DELAYTIME.
            elist1(list,CURTIME) = DELAYTIME;
        else
            DELAYTIME = getcurtime(node, elist2) + delay;
            list = find(elist2(:,CURTIME)-DELAYTIME < 0 & elist2(:,SRC)==node); 
            % Set the CURTIME field of all the rows in list to DELAYTIME.
            elist2(list,CURTIME) = DELAYTIME;
        end
            
	end
   
    % Move the first row of elist to SIMRESULT
    function updatesimlist()
        % update the simulation result after the transmission is complete
        % and remove the packet from the respective elist
        if bus1 
            SIMRESULT = [SIMRESULT; elist1(1,1:end)];
            % This command will delete the first row of elist1.
            elist1(1,:)=[];
        end
        
        if bus2
            SIMRESULT = [SIMRESULT; elist2(1,1:end)];
            % This command will delete the first row of elist2.
            elist2(1,:)=[];
        end
        
    end

    function routingpackets(src, dst)
        
        % Routers are referred with the index as R1 = 1; R2 = 2; R3 = 3; R4 = 4;
        % networknodes = [R1, R2, R3, R4];
        adj = [0 0 1 1; 0 0 1 1; 1 1 0 0; 1 1 0 0];               % adjacency list of routers in network; index corresponds to the resp router
        c1 = randi([1,10],1); c2 = randi([1,10],1);               % uniformly distribute cost of each link between [1, 10]
        c3 = randi([1,10],1); c4 = randi([1,10],1);               % c1 = w(r1, r3); c2 = w(r1, r4); c3 = w(r2, r3); c4 = w(r2, r4)
        edgeweights = [0 0 c1 c2; 0 0 c3 c4; c1 c2 0 0; c3 c4 0 0];
        
        % Run Dijkstra's algorithm to calculate the least cost path
        [costs, paths] = dijkstra(adj, edgeweights);
        pathlength = cellfun('length', paths);                  % number of hops in the shortest path  
              
        if((src == 1 || src == 2) && (dst == 3 || dst == 4))    % if source is connected to bus1 and destination is connected to bus2
            rcosts = costs(1:2, 3:4);
            [~, index] = min(rcosts(:));                        % get the minimum cost path from the costs array returned by Dijkstras algo
            [srouter, drouter] = ind2sub(size(rcosts), index);  % find the source and destination routers corresponding to min cost path
            drouter = drouter + 2;                              % get destination router from the original costs array
            rtdelay = 8;                                        % transmission delay between the routers (1Gbps)
            % update the total delay depending on the route taken by the packet
            tdelay = tdelay + ((pathlength(srouter, drouter) - 1) * (rtdelay + pd)) + (td + pd);         
        end
        
        if((src == 3 || src == 4) && (dst == 1 || dst == 2))    % if source is connected to bus2 and destination is connected to bus1
            rcosts = costs(3:4, 1:2);
            [~, index] = min(rcosts(:));                        % get the minimum cost path from the costs array returned by Dijkstras algo
            [srouter, drouter] = ind2sub(size(rcosts), index);  % find the source and destination routers corresponding to min cost path
            srouter = srouter + 2;
            tdelay = tdelay + ((pathlength(srouter, drouter)) * (td + pd));         
        end
    end

    function backoffoncollision
        % update the elist1 for collisions in bus1
        if(timediff1 <= pd)
            if elist1(1,TXTIME) == 0
                elist1(1,TXTIME) = elist1(1,CURTIME);
            end
            if elist1(2,TXTIME) == 0
                elist1(2,TXTIME) = elist1(2,CURTIME);
            end
            elist1(1, COLLISIONS) = elist1(1, COLLISIONS) + 1;
            elist1(2, COLLISIONS) = elist1(2, COLLISIONS) + 1;

            if (elist1(1, COLLISIONS) < maxbackoff)
                bk(src)=(randi(2^(elist1(1,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(src)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
            if((src == 1 && dst == 2) || (dst == 1 && src == 2))
                if (elist1(2,COLLISIONS)<maxbackoff)
                    bk(dst)=(randi(2^(elist1(2,COLLISIONS)),1,1)-1)*tbackoff;
                else
                    bk(dst)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
                end
                % Delay the packets at destination node only if it is connected to the same bus.
                delaynodepkts(dst, pd - timediff1 + bk(dst));
            end
            % Delay the packets at source node.
            delaynodepkts(src, pd + timediff1 + bk(src));
        end
        if(timediff2 <= pd)
            % update the elist2 for collisions in bus2
            if elist2(1,TXTIME) == 0
                elist2(1,TXTIME) = elist2(1,CURTIME);
            end
            if elist2(2,TXTIME) == 0
                elist2(2,TXTIME) = elist2(2,CURTIME);
            end
            elist2(1, COLLISIONS) = elist2(1, COLLISIONS) + 1;
            elist2(2, COLLISIONS) = elist2(2, COLLISIONS) + 1;

            if (elist2(1, COLLISIONS) < maxbackoff)
                bk(src)=(randi(2^(elist2(1,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(src)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
            if((src == 3 && dst == 4) || (dst == 3 && src == 4))
                if (elist2(2,COLLISIONS)<maxbackoff)
                    bk(dst)=(randi(2^(elist2(2,COLLISIONS)),1,1)-1)*tbackoff;
                else
                    bk(dst)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
                end
                % Delay the packets at destination node only if it is connected to the same bus.
                delaynodepkts(dst, pd - timediff2 + bk(dst));
            end
            % Delay the packets at source node.
            delaynodepkts(src, pd + timediff2 + bk(src));
        end
    end
    
disp(toc);

end

function projectdemo2()
    clear all;
    tic
    
    % NODE Id's
    A=1; B=2; C=3; D=4;
    nodes = [A, B, C, D];
    
    % Elist is an (n by 7) matrix where n is the number of packets
    % waiting in the queue of nodes A and B. We define a packet to be a
    % row in the elist and is of the form
    % [SRC DEST GENTIME TXTIME RXTIME CURTIME COLLISIONS]
    elist = [];
    
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
    CLOCK = 0;
 
    TOTALSIM = 30*10^3; % Total simulation time
    lambda = .5;
    frameslot = 50; % frame slot time (usec) <-- changed from 500us
    td = 80;% transmission delay on BUS (usec)
    pd = 10; % propagation delay on BUS
    tdelay = td + pd; % total delay incurred during a pkt transmission
	tbackoff = frameslot; % time slot (usec) for backoff algorithm
    maxbackoff = 3; % maximum backoff time is 2^3 frame slot
	
    % the time at which the last packet was generated at node A, B, C and D resp.
    GENTIMECURSOR=[0 0 0 0];
        
    % create packets for all the nodes
    [x, y] = size(nodes);
     for n = 1:y
         createpacket(n);
     end
    
    if size(elist,1) == 0
        disp('No packets to simulate');
        return;
    end
    
    % Collect the statistics in this array.
    SIMRESULT = [];
        
    while(1)
        
        % Update the clock.
        updateclock();     
        % Find the source node of the packet at the first row in elist.
        src = elist(1,SRC);
		
        % Randomly generating destination nodes and adding to the elist.
        destnodes = nodes;
        destnodes(src) = [];
        randomindex = randi(length(destnodes));
        dst = destnodes(randomindex);
        elist(1,DEST) = dst;
        
        timediff = elist(2,CURTIME) - elist(1,CURTIME);
        
		if timediff > pd
			%% No collision case
            % Set the tx time. This should be the time when the packet
            % is transmitted for the first time.
            if elist(1,TXTIME) == 0
                elist(1,TXTIME) = elist(1,CURTIME);
            end
            % Set the rx time.
            elist(1,RXTIME) = elist(1,CURTIME) + tdelay;

            updatesimlist();
			
            createpacket(src);
			
			% add tdelay to CLOCK.
            delaypkts(tdelay);
        else
            % Collision and so Backoff
            if elist(1,TXTIME) == 0
                elist(1,TXTIME) = elist(1,CURTIME);
            end
            if elist(2,TXTIME) == 0
                elist(2,TXTIME) = elist(2,CURTIME);
            end
            elist(1,COLLISIONS) = elist(1,COLLISIONS)+1;
            elist(2,COLLISIONS)=elist(2,COLLISIONS)+1;
            
            if (elist(1,COLLISIONS)<maxbackoff)
                bk(src)=(randi(2^(elist(1,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(src)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
            if (elist(2,COLLISIONS)<maxbackoff)
                bk(dst)=(randi(2^(elist(2,COLLISIONS)),1,1)-1)*tbackoff;
            else
                bk(dst)=(randi(2^(maxbackoff),1,1)-1)*tbackoff;
            end
			
			% Delay the packets at each node.
            delaynodepkts(src,pd+timediff+bk(src));
            delaynodepkts(dst,pd-timediff+bk(dst));
        end
		
        if min(GENTIMECURSOR) > TOTALSIM
            disp('Completed!');
            calcstat();
            break;
        end
    end
       
	   function calcstat()
        nodename = ['A', 'B', 'C', 'D']; 
        figure(1);
        figure(2);
        plotcount = 0;
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
                    disp(var1);
                    
                    % number of packets sent from node i to node j
                    var2 = length(var1);
                    fprintf('\n number of packets is %d \n', var2);
                    
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
                    
                    % access delay from node i to node j
                    var4 = var1(:, RXTIME) - var1(:, TXTIME);
                    figure(2);
                    accessdelay = sum(var4);
                    plottitle = strcat('Access delay from node ', num2str(i), ' to ', num2str(j));
                    subplot(4,3, plotcount);
                    plot(1:length(var4),var4);
                    xlabel('Packet sequence #');
                    ylabel('Delay in \mu sec');
                    title(plottitle);
                    
                    % end to end throughput for all pair of nodes
                    meanendtoend = mean(var3);
                    
                    % Average end to end throughput
                    avgthroughput = ((1000*8)/meanendtoend)*10^6;  % bits/sec
               end
            end
            % total queuing delay at node A, B, C and D
            queuedelay = queuedelay - 3*tdelay;

        end
    end
	
    % The clock slips to the RXTIME i.e., add delay time to CLOCK.
    function delaypkts(delay)
        CLOCK=CLOCK+delay;
        % list will have the row number of elist whose CURTIME field value
        % is less than CLOCK. Remember CLOCK is now the RXTIME. 
        list = find(((elist(:,CURTIME)-CLOCK) < 0));
        % Set the CURTIME field of all the rows in list to CLOCK.
        elist(list,CURTIME)=CLOCK;
    end
   
    function updateclock()
        % SORTROWS(elist,CURTIME) sorts the rows of elist in ascending
        % order for the column CURTIME.
        elist=sortrows(elist,CURTIME);
        % Set the clock to the CURTIME of the packet in the first row of
        % the elist since this is the packet that contends first for the
        % channel. 
        CLOCK=elist(1,CURTIME);
    end
    
    function pkt = createpacket(nodeid)
        % Find the inter-arrival time.
        interarvtime = round(frameslot*exprnd(1/lambda,1,1));
        
        % Find the birth time.
        GENTIMECURSOR(nodeid) = GENTIMECURSOR(nodeid) + interarvtime;
        
        % Create the packet. Unknown fields are set to 0.
        % [SRC= nodeid GENTIME=birthtime TXTIME=0, RXTIME=0
        % CURTIME=birthtime COLLISIONS=0]
        pkt = [nodeid 0 GENTIMECURSOR(nodeid) 0 0 GENTIMECURSOR(nodeid) 0];
        
        % Enqueue to the event list. 
        elist = [elist; pkt];
    end
	
	% Return the CURTIME of node
	function t=getcurtime(node)
		idx=find(elist(:,SRC)==node,1,'first');
        t=elist(idx,CURTIME);
    end
	
    function delaynodepkts(node,delay)
        DELAYTIME = getcurtime(node) + delay;
		% idx will have the row number of elist whose CURTIME field value
        % is less than DELAYTIME.
        list=find(elist(:,CURTIME)-DELAYTIME < 0 & elist(:,SRC)==node); 
		% Set the CURTIME field of all the rows in list to DELAYTIME.
        elist(list,CURTIME)=DELAYTIME;
	end
   
    % Move the first row of elist to SIMRESULT
    function updatesimlist()
        % elist(1,1:end) means row 1 and all columns.
        % Enqueue to the the first row of event list to SIMRESULT. 
        SIMRESULT=[SIMRESULT; elist(1,1:end)];
        
        % This command will delete the first row of elist.
        elist(1,:)=[];
    end
 
    
disp(toc);

end

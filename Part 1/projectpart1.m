function projectpart1()
    clear all;
    tic
    
    % NODE Id's
    A=1;
    B=2;
    
    % elist as an (n by 6) matrix where n is the number of packets
    % waiting in the queue of nodes A and B. We define a packet to be a
    % row in the elist and is of the form
    % [SRC GENTIME TXTIME RXTIME CURTIME COLLISIONS]
    elist=[];
    
    % Columns in the event list
    SRC=1;      % Id of the source node
    GENTIME=2;  % Time at which the packet is created at the source node
    TXTIME=3;   % Time at which the packet is transmitted
    RXTIME=4;   % Time at which the packet is received at the dest node
    % The current simulation time w.r.t the packet. You can think of
    % CURTIME to be the time at which the simulator had last visited this
    % packet.
    CURTIME=5;  
    COLLISIONS=6; % Number of collisions incurred by 
                  % this packet due to channel contention
    
    % Reference:
    % http://en.wikipedia.org/wiki/Discrete_event_simulation#Clock 
    % Clock:
    % The simulation must keep track of the current simulation time,
    % whatever measurement units are suitable for the system being modeled. 
    CLOCK=0;
 
    TOTALSIM=30*10^3; % Total simulation time
    lambda = 0.5;
    frameslot = 50; % frame slot time (usec) 
    td = 80;% transmission delay on BUS (usec)
    pd = 10; % propagation delay on BUS
    tdelay = td + pd; % total delay incurred during a pkt transmission
	tbackoff = frameslot; % time slot (usec) for backoff algorithm
    maxbackoff = 3; % maximum backoff time is 2^3 frame slot
	
    % the time at which the last packet was generated at node A and B resp.
    GENTIMECURSOR = [0 0];
        
    % Create packet for node A. Function createpacket() is defined below.
    createpacket(A); 
    createpacket(B);
    
    if size(elist, 1) == 0
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
		dst = mod(src,2)+1;
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
			
			% add tdelay to CLOCK. Check delaypkts() for more details.
            delaypkts(tdelay);
        else
            % Collision and so Backoff
            if elist(1,TXTIME)==0
                elist(1,TXTIME)=elist(1,CURTIME);
            end
            if elist(2,TXTIME)==0
                elist(2,TXTIME)=elist(2,CURTIME);
            end
            elist(1,COLLISIONS)=elist(1,COLLISIONS)+1;
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
        AtoB=SIMRESULT(SIMRESULT(:,SRC)==A,:);
        BtoA=SIMRESULT(SIMRESULT(:,SRC)==B,:);
        
        % Total packets sent
        AtoBnum=length(AtoB);
        BtoAnum=length(BtoA);
        
        % Queue delay
        queuedelaya=AtoB(:,RXTIME)-AtoB(:,GENTIME);
        queuedelayb=BtoA(:,RXTIME)-BtoA(:,GENTIME);
        queuedelaya=queuedelaya-tdelay;
        queuedelayb=queuedelayb-tdelay;
        
        figure;
        subplot(3,2,1)
        plot(1:size(queuedelaya,1),queuedelaya(1:end))
        axis([0 AtoBnum 0 max(queuedelaya)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Queue delay at node A');
        
        subplot(3,2,2)
        plot(1:size(queuedelayb,1),queuedelayb(1:end))
        axis([0 BtoAnum 0 max(queuedelayb)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Queue delay at node B');   
        
        % Access delay
        accessdelaya=AtoB(:,RXTIME)-AtoB(:,TXTIME);
        accessdelayb=BtoA(:,RXTIME)-BtoA(:,TXTIME);
        
        subplot(3,2,3)
        plot(1:size(accessdelaya,1),accessdelaya(1:end))
        axis([0 AtoBnum 0 max(accessdelaya)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Access delay at node A');
        
        subplot(3,2,4)
        plot(1:size(accessdelayb,1),accessdelayb(1:end))
        axis([0 BtoAnum 0 max(accessdelayb)])
        xlabel('Packet sequence #');
        ylabel('Delay in \mu sec');
        title('Access delay at node B');
        
        % Frame interval
        subplot(3,2,5)
        frameintba=BtoA(2:end,GENTIME)-BtoA(1:end-1,GENTIME);
        plot(1:BtoAnum-1,frameintba/1000)
        axis([0 BtoAnum 0 max(frameintba)/1000])
        xlabel('Packet sequence #');
        ylabel('Frame interval in msec');
        title('Frame intervals at node B');
        
        subplot(3,2,6)
        frameintab=AtoB(2:end,GENTIME)-AtoB(1:end-1,GENTIME);
        plot(1:AtoBnum-1,frameintab/1000)
        axis([0 AtoBnum 0 max(frameintab)/1000])
        xlabel('Packet sequence #');
        ylabel('Frame interval in msec');
        title('Frame intervals at node A');
        
        
        % Histogram to verify if the distribution is exp
        figure;
        subplot(2,1,1)
        hist(frameintab,60);
        xlabel('frame intervals in \mu sec');
        ylabel('# of frames');
        subplot(2,1,2)
        hist(frameintba,60);
        xlabel('frame intervals in \mu sec');
        ylabel('# of frames');
        
        % mean end to end time
        meanendtoenda=mean(AtoB(:,RXTIME)-AtoB(:,GENTIME));
        
        meanendtoendb=mean(BtoA(:,RXTIME)-BtoA(:,GENTIME));
        
        % Average end to end throughput
        Avgtha=((1000*8)/meanendtoenda)*10^6;  % bits/sec
        Avgthb=((1000*8)/meanendtoendb)*10^6;  % bits/sec
        avgth2=mean(8000.*(AtoB(:,RXTIME)-AtoB(:,GENTIME)).^(-1));
        
        fprintf('Total packets sent from node A=%d\n',AtoBnum);
        fprintf('Total packets sent from node B=%d\n',BtoAnum);
        fprintf('Average frame interval at node A=%d\n',mean(frameintab));
        fprintf('Average frame interval at node B=%d\n',mean(frameintba));
        fprintf('Average access delay at node A=%d\n',mean(accessdelaya));
        fprintf('Average access delay at node B=%d\n',mean(accessdelayb));
        fprintf('Average queue delay at node A=%d\n',mean(queuedelaya));
        fprintf('Average queue delay at node B=%d\n',mean(queuedelayb));
        fprintf('Average end to end throughput from A to B=%d\n',Avgtha);
        fprintf('Average end to end throughput from B to A=%d\n',Avgthb);
        fprintf('Simulation end time=%d\n',CLOCK);
    end
	
    % The clock slips to the RXTIME i.e., add delay time to CLOCK.
    function delaypkts(delay)
        CLOCK=CLOCK+delay;
        % It might so happen that the new packet at the SRC node might have
        % been generated when the previous packet was in flight. This new
        % packet cannot be transmitted immediately and hence has to wait
        % till the previos packet has reached the destination. 
        
        % list will have the row number of elist whose CURTIME field value
        % is less than CLOCK. Remember CLOCK is now the RXTIME. 
        list=find(((elist(:,CURTIME)-CLOCK) < 0));
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
        GENTIMECURSOR(nodeid) = GENTIMECURSOR(nodeid)+interarvtime;
        
        % Create the packet. Unknown fields are set to 0.
        % [SRC= nodeid GENTIME=birthtime TXTIME=0, RXTIME=0
        % CURTIME=birthtime COLLISIONS=0]
        pkt = [nodeid GENTIMECURSOR(nodeid) 0 0 GENTIMECURSOR(nodeid) 0];
        
        % Enqueue to the event list. In matlab you can append a row x to a
        % matrix X by using the command X=[X; x]
        elist = [elist; pkt];
    end
	
	% Return the CURTIME of node
	function t=getcurtime(node)
		idx=find(elist(:,SRC)==node,1,'first');
        t=elist(idx,CURTIME);
    end
	
    function delaynodepkts(node,delay)
        DELAYTIME = getcurtime(node)+delay;
		% idx will have the row number of elist whose CURTIME field value
        % is less than DELAYTIME.
        list = find(elist(:,CURTIME) - DELAYTIME < 0 & elist(:,SRC) == node); 
		% Set the CURTIME field of all the rows in list to DELAYTIME.
        elist(list,CURTIME) = DELAYTIME;
	end
   
    % Move the first row of elist to SIMRESULT
    function updatesimlist()
        % Enqueue to the the first row of event list to SIMRESULT. 
        SIMRESULT=[SIMRESULT; elist(1,1:end)];
        
        % Delete the first row of elist.
        elist(1,:)=[];
    end
 
    
disp(toc);

end

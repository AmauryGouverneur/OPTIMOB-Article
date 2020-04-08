function [meas_GF,cost_GF,avgCostHist,minCostHist] = greedy_forward_algo(n_measurements,T,n_part,n_draw,y,meas_1_j)
%%Simulated Annealing algorithm

if nargin < 4 && nargin > 1
    n_part = 250; %number of particles in the particle filter
    n_draw = 100; %number of draws in the MC
end
if nargin < 5 
    online = false ;
    meas_1_j = 0;
    y = 0;
else
    online = true ;
end

visualizationFlag=0;       % 0 => don't visualize bit frequencies
                           % 1 => visualize bit frequencies

verboseFlag=0;             % 1 => display details of each generation
                           % 0 => run quietly
convergenceFlag=0;         % 1 => plot convergence curve
                           % 0 => does not


% pre-generate two ?repositories? of random binary digits from which the  
% the masks used in mutation and uniform crossover will be picked. 
% maskReposFactor determines the size of these repositories.


% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,n_measurements);
maxFitnessHist=zeros(1,n_measurements);


% To identify copies in population
elite = [] ; 
if online 
    candidates = (1:T);
else 
    candidates = (0:T);
end

pop = candidates'; 
while length(elite)<n_measurements  
    
    
    % evaluate the fitness of the population. The vector of fitness values 
    % returned  must be of dimensions 1 x popSize.
    gen = length(elite)+1;
    popSize = length(candidates);
    fitnessVals=localFitnessFunction(pop);
    [maxFitnessHist(1,gen),maxIndex]=max(fitnessVals);
    avgFitnessHist(1,gen)=mean(fitnessVals,'omitnan');

    elite = pop(maxIndex,:);
     
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    % Conditionally perform bit-frequency visualization
    if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen),'%3.3f') ]);
    end
    if visualizationFlag
        figure(1)
        set (gcf, 'color', 'w');
        hold off
        if online
            histogram(pop,1:T,'Normalization','countdensity'); hold on;
            plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
            axis([1 T 0 popSize]);
        else 
            histogram(pop,0:T,'Normalization','countdensity'); hold on;
            plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
            axis([0 T 0 popSize]);
        end
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen))]);
        ylabel('Frequency of measure in t');
        xlabel('time t');
        drawnow;
    end    
    
    pop = [ones(length(candidates),length(elite)).*elite,candidates'];
    pop = sort(pop,2);
end

avgCostHist = -avgFitnessHist;
minCostHist = -maxFitnessHist;

meas_GF = elite;
cost_GF = -maxFitnessHist(end);

%% plot and print
if convergenceFlag
    figure
    set(gcf,'Color','w');
    hold off
    plot(1:n_measurements,-avgFitnessHist,'k-');
    hold on
    plot(1:n_measurements,-maxFitnessHist,'c-');
    title('Minimum and Average Cost');
    xlabel('Generation');
    ylabel('Cost');
end

    function fitness = localFitnessFunction(pop)
       % function to MAXIMIZE
       fitness = zeros(size(pop,1),1);
        %parfor j = 1:size(pop,1)
        for j = 1:size(pop,1)
            meas = pop(j,:);
            if online
                fitness(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part,y,meas_1_j);
            else 
                fitness(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part);
            end
        end
    end
end

function [meas_GF,cost_GF,avgCostHist,minCostHist] = greedy_forward_with_release_algo(n_measurements,T,n_part,n_draw,y,meas_1_j)
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

verboseFlag = 1;

convergenceFlag=1;         % 1 => plot convergence curve
                           % 0 => does not


% pre-generate two ?repositories? of random binary digits from which the  
% the masks used in mutation and uniform crossover will be picked. 
% maskReposFactor determines the size of these repositories.


% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,2*n_measurements);
maxFitnessHist=zeros(1,2*n_measurements);


elite = [] ; 
if online 
    initial_candidates = (1:T);
else 
    initial_candidates = (0:T);
end
candidates = initial_candidates;
pop = candidates'; 
while length(elite)<n_measurements  
    
    
    % evaluate the fitness of the population. The vector of fitness values 
    % returned  must be of dimensions 1 x popSize.
    gen = length(elite)+1;
    fitnessVals=localFitnessFunction(pop);
    [maxFitnessHist(1,gen),maxIndex]=max(fitnessVals);
    avgFitnessHist(1,gen)=mean(fitnessVals,'omitnan');

    elite = pop(maxIndex,:);
    candidates_measurement_times = ones(1,length(initial_candidates)); 
    candidates_measurement_times(elite+1) = 0;
    candidates = initial_candidates(ones(1,length(initial_candidates)) == candidates_measurement_times);
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    % Conditionally perform bit-frequency visualization
    
     if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen),'%3.3f') ]);
    end
    pop = [ones(length(candidates),length(elite)).*elite,candidates'];
    pop = sort(pop,2);
end

maxFitnessHist(1,gen) = MC_MSE_estimator(elite,T,1000,1000);
release_max = n_measurements;
release = 1;

while release <= release_max
    picked_out = ceil(n_measurements*rand());
    elite_released = elite(cat(2, 1:picked_out-1, picked_out+1:n_measurements));
    
    candidates_measurement_times = ones(1,length(initial_candidates)); 
    candidates_measurement_times(elite_released+1) = 0;
    candidates_released = initial_candidates(ones(1,length(initial_candidates)) == candidates_measurement_times);
    
    pop = [ones(length(candidates_released),length(elite_released)).*elite_released,candidates_released'];
    pop = sort(pop,2);
    
    fitnessVals=localFitnessFunction(pop);
   
    [maxFitnessHist(1,gen+release),maxIndex]=max(fitnessVals);
    avgFitnessHist(1,gen+release)=mean(fitnessVals,'omitnan');
    
    fitness_new_elite = MC_MSE_estimator(pop(maxIndex,:),T,1000,1000);
    if(fitness_new_elite  < maxFitnessHist(1,gen+release-1))
        maxFitnessHist(1,gen+release) = maxFitnessHist(1,gen+release-1);
    else
        maxFitnessHist(1,gen+release) = fitness_new_elite;
        elite = pop(maxIndex,:);
    end
    
    if verboseFlag
        display(['gen=' num2str(gen+release,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen+release),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen+release),'%3.3f') ]);
    end
    
    
    release = release + 1 ; 
end




avgCostHist = -avgFitnessHist;
minCostHist = -maxFitnessHist;

meas_GF = elite;
cost_GF = MC_MSE_estimator(elite,T,1000,1000);

%% plot and print
if convergenceFlag
    figure
    set(gcf,'Color','w');
    hold off
    plot(1:length(avgFitnessHist),-avgFitnessHist,'k-');
    hold on
    plot(1:length(avgFitnessHist),-maxFitnessHist,'c-');
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

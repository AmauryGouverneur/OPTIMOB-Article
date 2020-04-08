function [meas_SA,cost_SA,avgCostHist,minCostHist] = SA_algo(n_measurements,T,pop_size,max_gen,n_part,n_draw,y,meas_1_j)
%%Simulated Annealing algorithm

if nargin < 6 && nargin > 4
    n_part = 250; %number of particles in the particle filter
    n_draw = 100; %number of draws in the MC
end
if nargin < 8 
    online = false ;
    meas_1_j = 0;
    y = 0;
else
    online = true ;
end

len=n_measurements;        % The length of the genomes  
popSize=pop_size;          % The size of the population (must be an even number)
maxGens=max_gen;                % The maximum number of generations allowed in a run
probMutation=0.3;        % The mutation probability (per bit)

visualizationFlag=0;       % 0 => don't visualize bit frequencies
                           % 1 => visualize bit frequencies

verboseFlag=0;             % 1 => display details of each generation
                           % 0 => run quietly
convergenceFlag=0;         % 1 => plot convergence curve
                           % 0 => does not

useMaskRepositoriesFlag=0; % 1 => draw uniform mutation masks from 
                           %      a pregenerated repository of randomly generated bits. 
                           %      Significantly improves the speed of the code with
                           %      no apparent changes in the behavior of
                           %      the SSA
                           % 0 => generate uniform crossover and mutation
                           %      masks on the fly. Slower.

% pre-generate two ?repositories? of random binary digits from which the  
% the masks used in mutation and uniform crossover will be picked. 
% maskReposFactor determines the size of these repositories.

maskReposFactor=5;
mutmaskRepos=rand(popSize,(len+1)*maskReposFactor)<probMutation;
temperature = 20;

% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,maxGens+1);
maxFitnessHist=zeros(1,maxGens+1);
rateAcceptanceHist=zeros(1,maxGens+1);

% the population is a popSize by len matrix of randomly generated boolean
% values

pop = zeros(popSize,len);
for i=1:popSize
    if online 
        pop(i,:) = sort(randperm(T,len));
    else 
        pop(i,:) = sort(randperm(T+1,len)-1);
    end
end

% To identify copies in population
pop = sortrows(pop);
gen = 0 ; 

while gen<=maxGens  
    % evaluate the fitness of the population. The vector of fitness values 
    % returned  must be of dimensions 1 x popSize.
    counter = 0; 
    fitnessVals=localFitnessFunction(pop);
    
    
    if gen>0 && maxFitnessHist(1,gen) > max(fitnessVals)
        maxFitnessHist(1,gen+1)= maxFitnessHist(1,gen);
    else 
        [maxFitnessHist(1,gen+1),maxIndex]=max(fitnessVals);
    end
    avgFitnessHist(1,gen+1)=mean(fitnessVals,'omitnan');
     
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    % Conditionally perform bit-frequency visualization
    
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
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen+1))]);
        ylabel('Frequency of measure in t');
        xlabel('time t');
        drawnow;
    end
    
    % implement mutations
    if useMaskRepositoriesFlag
        temp=floor(rand*len*(maskReposFactor-1));
        masks=mutmaskRepos(:,temp+1:temp+len);
    else
        masks=rand(popSize, len)< probMutation*ones(popSize,len);
    end
    % masks(i,j)==1 iff pop(i,j) has to be mutated (0 elsewhere)
    if online 
        pop_mutation = sort((1-masks).*pop + masks.*(unidrnd(T,popSize,len)),2);
    else 
        pop_mutation = sort((1-masks).*pop + masks.*(unidrnd(T+1,popSize,len)-1),2);
    end
    % Replace duplicates measurements and sort
    %parfor i=1:popSize
    for i=1:popSize
        if online 
            pop_mutation(i,:) = replace_duplicates(pop_mutation(i,:),T);
        else 
            pop_mutation(i,:) = replace_duplicates(pop_mutation(i,:),T);
        end
    end
    pop_mutation = sort(pop_mutation,2);
    pop_mutation = sortrows(pop_mutation);
    
    fitnessVals_mutation =localFitnessFunction(pop_mutation);
    
    for i=1:popSize
        if fitnessVals_mutation(i)>fitnessVals(i)
            pop(i,:) = pop_mutation(i,:);
            counter = counter+1;
        else
            if rand < exp(-(fitnessVals(i)-fitnessVals_mutation(i))/temperature)
                if (pop(i,:)~=pop_mutation(i,:))
                    pop(i,:) = pop_mutation(i,:);
                    counter = counter+1;
                end
            end
        end
    end
    
    
    pop = sort(pop,2);
    pop = sortrows(pop); % to identify copies in population
    
    temperature = temperature*0.8;
    rateAcceptanceHist(gen+1)=counter/pop_size;
    if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen+1),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen+1),'%3.3f') '   rateAcceptanceHist=' ...
            num2str(rateAcceptanceHist(1,gen+1),'%3.3f') '   temperature=' num2str(temperature,'%3.3f')]);
    end
    gen = gen+1; 
end

display(['fraction of changes= ' num2str(sum(rateAcceptanceHist)/maxGens,'%3.3f')]);
avgCostHist = -avgFitnessHist;
minCostHist = -maxFitnessHist;

meas_SA = sort(pop(maxIndex,:));
cost_SA = -maxFitnessHist(end);


%% plot and print
if convergenceFlag
    figure
    set(gcf,'Color','w');
    hold off
    plot(0:maxGens,-avgFitnessHist,'k-');
    hold on
    plot(0:maxGens,-maxFitnessHist,'c-');
    title('Minimum and Average Cost');
    xlabel('Generation');
    ylabel('Cost');
end

    function fitness = localFitnessFunction(pop)
       % function to MAXIMIZE
       % is designed to compute only once the fitness in cases of copies of
       % individuals
        [popLocSize,~] = size(pop);
        
        % first time that an individual appears
        indFirstCopy = find(sum( (pop(1:end-1,:)-pop(2:end,:)).^2 ,2)~=0)'+1;
        indFirstCopy = [1 indFirstCopy];
        
        % measurements of these first individuals
        firstMeas = pop(indFirstCopy,:);
        
        firstFitnesses = zeros(1,length(indFirstCopy)-1);
        %parfor j = 1:length(indFirstCopy)
        for j = 1:length(indFirstCopy)
            meas = firstMeas(j,:);
            if online
                firstFitnesses(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part,y,meas_1_j);
            else 
                firstFitnesses(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part);
            end
        end
        
        % copy the fitnesses for similar individuals
        fitness = repelem(firstFitnesses,diff([indFirstCopy popLocSize+1]));
    end
end

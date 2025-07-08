function [t,x,y,z,M,F,N,wm,wf,resistance,tau] = TargetSiteResistanceSims_while(K,N0,psi,x0,y0,Wm,Wf,Mu,stoichiometry,driveheterozgotes,Rm,T,ratio,Gauss, test_allele, test_threshold, Cswitch)


% MIT License
% 
% Copyright (c) 2021 Bhavin S Khatri
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Function that runs simulations given various input parameters. The degree
% of multiplex is implicit in the number of alleles defined and the
% specification of the fitness and mutation matrices, as well as the
% stochiometry and driveheterozygotes arrays

% K - number of alleles 
% Wm - male fitness matrix
% Wf - females fitness matrix
% N0 - initial population size
% Mu - matrix of mutation rates between kth and jth allele - size(Mu) = KxK
% stochiometry — vector of the change in fractions of alleles due to
%   action of drive
% driveheterozygotes — 2xNhz matrix whose columns are the allele index
%   corresponding to those genotypes involved in conversion of wild type ->
%   drive in gametogenesis
% T - number of generations


%Various options for Gauss
%if Gauss = 0 — multinomial random numbers — together with Cswitch=1 will use
%GSL version of multinomial random numbers, which is much quicker — but needs
%compilation on given platform
%if Gauss = 1 — Gaussian-Poisson approx for multinomial, where for alleles
%below threshold are treated with Poisson, and above Gaussian


%Threshold for the mean frequency of any of alleles  - if any of the
%alleles is below this threshold then the GSL version of multinomial is
%called
Gamma = 10;
%Turn population dynamics on or off
BH=1;

% err=0;

if numel(x0)~=K | numel(y0) ~= K
    disp('Length of initial vector x0,y0 should be = number of alleles K')
%     err=1;
    return
end

%Each column of stoichiometry S is vector describing the relative change in
%frequencies of each allele due to non-Mendelian drive
S = stoichiometry;

% Each column of driveheterozgotes H is a 2-element vector   
% containing the allele types involved in the heterozygotes involved in drive 
% (e.g. W\D for single target site, WW\DD, WR\DD,WN\DD for duplex target sites
% so H(:,1) = [1;2] represents allele W and allele D for single target site) 
H = driveheterozgotes;

%Above the number of columns represents the total number of heterozygotes
%involved in drive "reactions"
[~,Nhz] = size(H);

%Density dependent parameter from Beverton-Holt model of population
%dynamics

% mean_constraint = 'geometric'

phi=0;
omega = 2*pi/12; %This assumes that each generation is 1 month, so the period of oscillations is 12 months
% ratio = 100; %N1/N2
%If ratio isn't empty then assume N0 is the geometric mean of max and min
%i.e. N0 = sqrt(N1*N2) => N0 = sqrt(ratio*N2^2) = N2*sqrt(ratio) => N2 =
%N0/sqrt(ratio)

%For an arithemetic constaint N0=(N1+N2)/2 = N1*(1+ratio)/2 => N1 =
%2*N0/(1+ratio)
% ratio
if ~isempty(ratio)
    
%     N2 = N0/sqrt(ratio);
%     N1 = ratio*N2;

    N2 = 2*N0/(1+ratio);
    N1 = ratio*N2;

else
    alpha = N0/(Rm-1);
end

%Intial numbers of males and females
M0 = psi*N0;
F0 = (1-psi)*N0;


%Intialise population vectors
N = zeros(1,T+1);
M = zeros(1,T+1);
F = zeros(1,T+1);

N(1,1) = N0;
M(1,1) = M0;
F(1,1) = F0;

%intial frequency vectors of males and females - note these are not on
%simplex - so sum at any time point should equal 1


z0 = [x0*M0/N0;y0*F0/N0];

% size(x0)
% size(y0)
% size(z0)

%Intialise frequency vectors
x = zeros(K,T+1);
y = zeros(K,T+1);
z = zeros(2*K,T+1);

x(:,1) = x0;
y(:,1) = y0;
z(:,1) = z0;

%Intialise mean fitness vectors
wwf = zeros(1,T+1);
wwm = zeros(1,T+1);


% N0 = alpha*(Rm-1);


U=0;
k=2;

tau=NaN;


while U~=1
        
%         tic
%         %Generation time t =k-1
%         disp(['t = ',num2str(k-1)]);



%% Mean fitness in generation t
        %frequency vectors in generation t
        xt = x(:,k-1);
        yt = y(:,k-1);
        
         
        %Mean fitness of males and females in generation t
        wm = yt'*Wm*xt;
        wf = xt'*Wf*yt;
        

        wwm(k-1) = wm;
        wwf(k-1) = wf;

        
        
        
        
%% Calculate new expected population sizes in generation t+1 using modified Beverton-Holt model
        %Population sizes in generation t
%         disp(['t=',num2str(k-1)]);
        Mt = M(k-1);
        Ft = F(k-1);
        Nt = N(k-1);
        
        % alpha

        if BH==0
            NN=Nt;
            FF=NN/2;
            MM=NN/2;
        else

            if ~isempty(ratio)
                Nosc = (N1+N2)/2 + (N1-N2)/2*sin(omega*(k-1) + phi)
                alpha = Nosc./(Rm-1);
            end

            if Ft==0
                FF = 0;
                MM = 0;
            else
                FF = Rm*wf*Ft/(1+Nt/alpha);
    %             MM = Rm*wm*Ft/(1+Nt/alpha);
                MM = Rm*wf*Ft/(1+Nt/alpha); 
            end
            
            
            
            %This is mean expected population size — draw below from Poisson to get
            %realised population size
            NN = MM+FF;
    
            
    
            
            %% Draw a new random NN from Poisson disitribution
           NN = poissrnd(NN); %Must be an integer
        

       end


       if NN~=0
           
%% Mean/expected frequencies in generation t+1



        %Update allele frequencies 
        [xx, yy] = AlleleFrequencyUpdate(xt,yt,Wm,Wf,wm,wf,Mu,H,S,K);



        %Population level frequencies of male and female alleles
        zz = [xx/2; yy/2];
                
                
        if ~isempty(find(zz<0, 1))
            disp('Error: frequency vector must not have negative entries')
            zz
            xx
            yy
            xt
            yt
            x0
            y0
            k
            NN
            N0
            Gauss
            'hello'
            N(k-10:k-1)


        end



        nn = WrightFisherSampling(zz,NN,Gauss,Cswitch);

        %Population level frequencies in generation t+1
        z(:,k) = nn/NN;


        if ~isempty(find(z(:,k)<0, 1))
            disp('Error: frequency vector must not have negative entries')
            save('temp.mat')
            z(:,k)
            nn
            zz
            xx
            yy
            xt
            yt
            x0
            y0
            k
            NN
            N0
            Gauss
             'hello'
            N(k-10:k-1)
           

        end


        if BH==0
            M(k) = MM;
            F(k) = FF;
            N(k) = NN;
        else

        %Male & female populations in generation t+1
        M(k) = sum(nn(1:K));
        F(k) = sum(nn(K+1:2*K));
        N(k) = NN;

        end

        %Frequencies of alleles in amongst males and females in generation
        %t+1

        if M(k) ==0
            x(:,k) = zeros(K,1);
        else
            x(:,k) = nn(1:K)/M(k);
        end

        if F(k) ==0
            y(:,k)=zeros(K,1);
        else
            y(:,k) = nn(K+1:2*K)/F(k);
        end
            
    
            

    
%% Test end condition 
% (U=0 represents normal state; U=1 ends the simulations; 
%  U=2 indicates resistance frequency threshold has been reached, 
%  record this and let simulations run for 50% longer
%  N.B. Number of generations is k-1

        xr = sum(x(test_allele,k));
        yr = sum(y(test_allele,k));
%         U

        if xr >= test_threshold & yr>= test_threshold & U~=2


            U=2;
            T = round(k*1.5); %even though resistance condition has been met, run simulations for 50% longer to see behaviour
            resistance = 1;
            tau=k-1;

            k = k+1;

        elseif k==T+1

             if U~=2
                resistance = 0;
             end      

            U=1;

        else
            k=k+1;
        end
            
            
            

%         if numel(test_allele)<2
% 
%                 if x(test_allele,k) >=test_threshold & y(test_allele,k)>=test_threshold & U~=2
% 
% 
%                     U=2;
%                     T = round(k*1.5);
%                     establishment = 1;
% 
% 
%                     k=k+1;
% 
%                 elseif k==T+1
% 
%                     if U~=2
%                         establishment = 0;
%                     end
% 
%                     U=1;
% 
%                 else
%                     k=k+1;
%                 end
% 
%         elseif numel(test_allele)==2
%             
%                 if (x(test_allele(1),k) >=test_threshold & y(test_allele(1),k)>=test_threshold) & U~=2 | (x(test_allele(2),k) >=test_threshold & y(test_allele(2),k)>=test_threshold) & U~=2
% 
%                     U=2;
%                     T = round(k*1.5);
%                     establishment = 1;
% 
% 
%                     k=k+1;
% 
%                 elseif k==T+1
% 
%                     if U~=2
%                         establishment = 0;
%                     end
% 
%                     U=1;
% 
%                 else
%                     k=k+1;
%                 end
% 
% 
%         elseif numel(test_allele)>2
% 
%             xr = sum(x(test_allele,k));
%             yr = sum(y(test_allele,k));
% 
%             if xr >test_threshold & yr>test_threshold & U~=2
% 
% 
%                 U=2;
%                 T = round(k*1.5);
%                 establishment = 1;
% 
%                 k = k+1;
% 
%             elseif k==T+1
% 
%                  if U~=2
%                     establishment = 0;
%                  end      
% 
%                 U=1;
% 
%             else
%                 k=k+1;
%             end
% 
%         end
            
        


        else
            U=1;
            resistance = 0;
            
        end
%         toc
end


% 
t = 0:k-1;%T;
% T
x = x(:,1:k);
y = y(:,1:k);
z = z(:,1:k);

M = M(1:k);
F = F(1:k);
N = N(1:k);

% k
% size(t)
% size(M)
% size(wwm)
% size(wwf)

x(:,M==0) = NaN;
y(:,F==0) = NaN;



% wm = wwm(1:k-1);
% wf = wwf(1:k-1);
if numel(wwm)<k
    wm = wwm(1:k-1);
    wm = [wm,NaN];
else
    wm = wwm(1:k);
end

if numel(wwf)<k
    wf = wwf(1:k-1);
    wf = [wf,NaN];
else
    wf = wwf(1:k);
end


% size(wm)
% size(wf)
wm(M==0 |F==0) = NaN;
wf(F==0|M==0) = NaN;

% size(wm)
% size(wf)
% 
% 'end'
   

end


function [xx, yy] = AlleleFrequencyUpdate(xt,yt,Wm,Wf,wm,wf,Mu,H,S,K)


% MIT License
% 
% Copyright (c) 2021 Bhavin S Khatri
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% Update allele frequencies in next generation using scheme specified in
% scheme. 
% Selection is viability selection on zygote to adult stage 
% Mutation corresponds to spontaneous de novo mutation introduced due to
% error in replication. This is contained in matrix Mu whose elements are
% **** Mu(i,j) = rate(j->i) ****
% Drive is homing and an errors introduced due to drive action like errors
% in NHEJ

%N.B. these update schemes assumes mutation is weak (mu<<1 - if not then drive terms will have mutation dependence)
%Above the number of columns represents the total number of heterozygotes
%involved in drive "reactions"
[~,Nhz] = size(H);


    
       if wm==0
           xx = zeros(K,1);
       else
           xx = 1/2*(eye(K)+Mu)*((Wm*xt).*yt +(Wm*yt).*xt)/wm;
       end
       
       if wf==0
           yy = zeros(K,1);
       else
           yy = 1/2*(eye(K)+Mu)*((Wf*yt).*xt +(Wf*xt).*yt)/wf;
       end
       
        
        
        drivem = zeros(K,1);
        drivef = zeros(K,1);
        
        for j=1:Nhz
            
            %index of alleles involved in each heterozygote
            h1 = H(1,j);
            h2 = H(2,j);
            
            %Frequency of that heterozygote in generation t
            %**N.B. this is not the frequency of the heterozygote
            %calculated after selection, i.e. xw(h1)*yw(h2)+xw(h2)*yw(h1), 
            %assuming random fusion of gametes
            %since the number of conversions due to drive is proportional
            %the number of heterozygotes after selection 
            %[The W(h1,h2)/w term below accounts for selection in this heterozygote], 
            %but before gametogenesis
%             xt(h1)
%             yt(h1)
%             xt(h2)
%             yt(h2)
            
            hzfreq = xt(h1)*yt(h2) + xt(h2)*yt(h1);
            mut  = 1+Mu(h1,h1); %This only works if there is no backmutation back into this heterozygote, 
            %such that there will be drive terms dependent on the frequency
            %of other genotypes. In other words the drive term is always
            %dependent on.
            %N.B. diagonal terms of Mu are always negative as they are net
            %rate out of state, so above is basically 1- mu
            
            if wm==0
                drivem = drivem +zeros(K,1);
            else
                drivem  = drivem + 1/2*mut*S(:,j)*Wm(h1,h2)/wm*hzfreq;
            end
            
            if wf==0
                drivef = drivef + zeros(K,1);
            else
                drivef  = drivef + 1/2*mut*S(:,j)*Wf(h1,h2)/wf*hzfreq;
            end
            
            
        end
        
        xx = xx + drivem;
        yy = yy + drivef;
        
        
        
end




function nn = WrightFisherSampling(zz,NN,Gauss,Cswitch)



    %% Calculate new random populations and frequencies in t+1 accounting for drift

            % Draw NN individuals with frequency vector zz using
            % multinomial/Gaussian distribution (in vector nn) - mnrnd always
            % produces a row vector

%             zz=zz/sum(zz);

            if ~isempty(find(zz<0, 1))
                disp('Error: frequency vector must not have negative entries')
                zz
                
            else
                zz=zz/sum(zz);

                
                if abs(sum(zz)-1)>3*eps
%                     k
                    zz
                    norm(zz,1)-1
                    disp('Error: elements of frequency vector must sum to 1')
                end
            
            end
    %         sum(zz)-1


            if Gauss ==0
                
                if Cswitch==1
                    nn = Multinomial_gsl(NN,zz,numel(zz));
                    nn = double(nn);
                else
                    nn = mnrnd(NN,zz);
                end
                
                nn = nn';
                
                
            else
                                
                nn = GaussPoissonHybrid_mnrnd(NN,zz);

            end



            
            
            
end



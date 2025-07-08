function [resistance,Nend,tau,N,M,F,x,y,t] = TargetSiteResistance(m,s,psi,hDW,hDP,hN,sigmaR,sigmaP,Rm,K,xinit,yinit,epsilonW,epsilonP,nuW,nuP,muSite,betaWR,betaWP,betaPR,betaPW, xiR,xiP,preexisting,T,test_threshold,Gauss,Cswitch,plotfig,plotmale,vis,name)


% MIT License
% 
% Copyright (c) 2025 Bhavin S Khatri
% % 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% % furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.



% This code runs a single instance of Wright-Fisher pop-gen and Beverton-Holt population dynamics simulations of 
% a diecious population for m gRNAs, where each cleavage target site has 5 alleles: W - wild type; 
% R (R1) - functional resistant; N (R2) — non-functional
% resistant; P (R3) - partial resistant allele; D - drive, but where a mosaic of drive on a single chromosome
% is not possible (e.g. gametes WRD is not possible, but only DDD)

%Inputs:
% m — number of gRNAs (can in principle do simulations for any number, but
% computation becomes slower)
% psi — initial fraction of males in population 
% s — female fitness cost of homozygote with deleterious haplotypes on both chromosomes
% hDW — dominance coeff/fitness cost of each wild type allele in haplotype
% with drive on other chromosome
% hDP — dominance coeff/fitness cost of each partial resistance (R3) allele in haplotype
% with drive on other chromosome
% hN — dominance coeff/fitness cost of each non-functional resistance allele (R2) in haplotype 
% with drive on other chromosome
% sigmaR — female fitness cost of each occurrence of R (R1) in a genotype (e.g.
%   w(WWR/WWR) = (1-σR)^2, w(WWW/WWR) = 1-σR, w(WRR/RRR) = (1-σR)^5 )
% sigmaP — female fitness cost of each occurrence of P (R3) in a genotype (e.g.
%   w(WWP/WWP) = (1-σP)^2, w(WWW/WWP) = 1-σP, w(WPP/PPP) = (1-σP)^5 )
% Rm — absolute or intrinsic growth rate of population assuming population
%   is all WWW
% K — carrying capacity of population assuming population is all WW...W
% xinit — initial frequency vector of alleles in males
%         if scalar then assumes it is =xD — initial frequency of drive in males (DD...D)
% yinit — initial frequency vector of alleles in females 
%         if scalar then assumes it is =yD — initial frequency of drive in females (DD...D)

% epsilonW — efficiency of drive cleavage per target site with W allele
% epsilonP — efficiency of drive cleavage per target site with P allele
% nuW — non-homologous end-joining (NHEJ) rate per generation per individual
%   per target site with W allele
% nuP — non-homologous end-joining (NHEJ) rate per generation per individual
%   per target site with P allele
% muSite — mutation rate per cleavage target site per generation per individual
%   per target site (e.g. if a target site contains 10bp and bp mutation rate is 1e-9, set mu=1e-8)
% betaXY — fraction of type Y alleles from NHEJ mutations from allele X
% (i.e.. betaWR are the fraction of NHEJ mutations from W that give R
%        betaWP are the fraction of NHEJ mutations from W that give P,
%        betaPR are the fraction of NHEJ mutations from P that give R,
%        betaPW are the fraction of NHEJ mutations from P that give W

% xiR — fraction of functional resistant single nucleotide mutations that
% produce R (R1) alleles from W
% xiP — fraction of functional resistant single nucleotide mutations that
% produce P (R3) alleles from W

% preexisting — preexisting =1 runs simulations for a period 1/sigma as a
%   burn-in period to establish a mutation-selection balance
% T — length of simulations in generations (simulations will end earlier if population eliminated)
% test_threshold — at what frequency should the sum of all resistance
% alleles reach to terminate simulations (and signify resistance=1)
% Gauss — Gauss = 0 uses multinomial distribution; Gauss = 1 uses
%   Poisson-Gaussian hybrid approximation to the multinomial distribution.
%   If Gauss=0 the maximum population size is ~N=1e9 and is very slow — if
%   Gauss=1, the hybrid approx can be used with arbitrarily large N with
%   little speed penalty
% Cswitch — Cswitch = 1 if together with Gauss=0 uses GSL C version of
%   multinomial random number generator, which is much quicker than the
%   matlab implementation. Cswitch = 0 uses Matlab's multinomial RNG,
%   unless Gauss=1
% plotfig — plotfig=1 plots allele frequency time series and population
%   dynamics; plotfig=0 doesn't
% plotmale — plots male allele frequencies or not
% vis — vis=0 with plotfig=1 plots figures and makes them invisible and 
%   saves them (for remote jobs); their visibility state will be 0 and so when
%   opened in Matlab this needs to be changed, e.g. openfig('Figure.fig','visible')
% name — name of directory to create and store figures 

%Outputs:
% resistance — resistance=0 or 1, whether simulation produced resistance 
%   (freq(R)>test_threshold); 
%   if resistance alleles have not passed threshold by T, then resistance=0
% Nend — final population size at the end of simulation
% tau — the time when resistance occured (freq(R)>test_threshold) or = NaN
%   if population eliminated
% N — total population size over time
% M — male population size over time
% F — female population size over time
% x — matrix whose columns are the frequency of each allele in males at subsequent
%   time points
% y — matrix whose columns are the frequency of each allele in females at subsequent
%   time points
%   For both of these the order of alleles is W R N P D — so the 3rd row of x 
%   is the time series of the freq(N)    
% t — corresponding vector of times (generations)

%N.B. note that in paper frequency of males have frequency y and females x

%For m=3 (triplex) need to decided how to order our super-alleles/haplotypes
%The ordering is chosen such that we take the  the outer/kronecker product 
%(W,R,N)*(W,R,N)*(W,R,N) 
%and then from the "upper diagonal" (i<=j<=k) of the resulting 3D array 
%we take the column order of these elements of this rank deficient matrix to give
%the super allele ordering:
%[WWW WWR WRR RRR WWN WRN RRN WNN RNN NNN]
% then including drive:
%[WWW WWR WRR RRR WWN WRN RRN WNN RNN NNN, DDD]



%Fraction of males
% psi=0.5;

% test_allele are RRR RRN, RNN and NNN, which are resistant — their indices
% are calculated below
% test_threshold = 0.95;

%Density dependent parameter from Beverton-Holt model of population
%dynamics
alpha = K/(Rm-1);

% demographic seasonality 
ratio=[];

% %Dominance coefficients
% hDW = h;

% e=1e-14;
% if s==1
%     s=1-e;
% end
% 
% hN = log(1-hN*s)/log(1-s);


if numel(sigmaR)==1
    sigmaR = sigmaR*ones(1,m);
elseif numel(sigmaR)~=m
    disp('Error: sigmaR should be a scalar or number of elements equal to m')
    return
end

if numel(sigmaP)==1
    sigmaP = sigmaP*ones(1,m);
elseif numel(sigmaP)~=m
    disp('Error: sigmaP should be a scalar or number of elements equal to m')
    return
end


% %Multiplex level
% m = 3;

%Number of alleles is n=5 — [W R N P D] 
% The R,N,P correspond to R1, R2, R3 resistance alleles
%site
n=5;

% %Number of super alleles (haplotypes) ns
% %This is the number of combinations (unordered) of 3-alleles  (W,R,N)
% %across m sites & is given by n^(m)/m! - where n^(m) is the rising
% %factorial
% % + 1 drive allele, since drive replaces all allelles at each site
% ns =  SumTensorDiagonal(n-1,m) +1;

%Number of ordered alleles across m target sites
na = (n-1)^m +1;

% Create logical linear index which correspond to at least a single
% instance of drive (allele n) across all combinations of ordered alleles.
indDrive = zeros(1,n^m);


for k=1:n^m
    u = linearind2sub(m,n,k);
    
    if numel(find(u==n))>0
        indDrive(k) = 1;
    else
        indDrive(k) = 0;
    end
    
%     U = [U;u];
end

indDrive = logical(indDrive);

%% Calculate genotype fitness matrices
Wm = ones(na);
% Wm(4,:) = 3*ones(1,na);
% Wm(:,4) = 3*ones(na,1);
% Wm
Wf = zeros(na);

% Heuristic dominance model



% W(i,j) will be a matrix of the fitness of the homozgote w(ij\ij) - where
% i and j specify haplotypes/super-alleles

W = zeros(na,1);
H = zeros(na,1);

for k=1:na-1
    
    %linear index is k
    
    %convert k to subscript indices contained in u
    u = linearind2sub(m,n-1,k);
    
    
    
    % nn is a (n-1)x1 vector of the number of times that allele appears in index u 
    for kk=1:n-1
        nn(kk) = numel(find(u==kk));
    end
    
    
    if nn(3) == 0 %i.e. there are no N alleles in haplotype u
        nW(k) = nn(1);
        nR(k) = nn(2);
        nP(k) = nn(4);
    else %i.e. even if there are some functional alleles W and R if there are any N 
         %then it is effectively deleterious — signify by NaN
        nR(k) = NaN;
        nP(k) = NaN;
    end
    
    w = 1;
    
%     if sum(nn([3,end])) == 0 %No deleterious
    
    if nn(3) == 0 %No deleterious
        
        % Calculate the "fitness" of the haplotype — i.e. notional fitness
        % if organism was haploid
        for kk=1:numel(u) %kk is looping over target sites
            if u(kk) == 2
                w=w*(1-sigmaR(kk));
            elseif u(kk) ==4
                w=w*(1-sigmaP(kk));
            end
        end
    
    
        
        
        W(k) = w;
%         H(k) = hR;
    else
        W(k) = 1-s; %If sum(nn(3:end))>0 then at least 1 N,D
        
    end
            
% %All this below should be after loop (but don't think it matters!)                
% W(na) = 1-s;
%     
% %Add last element refering to D^m haplotype 
% nR = [nR;NaN];
% nP = [nP;NaN];
% nW = [nW;NaN];
    
    
end

W(na) = 1-s;
    
[sx,sy]=size(nR);
%Add last element refering to D^m haplotype 

if sy==1
    nR = [nR;NaN];
    nP = [nP;NaN];
    nW = [nW;NaN];
else
    nR = [nR,NaN];
    nP = [nP,NaN];
    nW = [nW,NaN];
end



for i=1:na
    for j=1:na
%         na
%         [i,j]
%         
%         [nW(i),nR(i),nP(i)]
%         [nW(j),nR(j),nP(j)]
        
        if isnan(nR(i)) && isnan(nR(j)) %is this not redundant??
%             disp('drive on i and j')
            Wf(i,j) = 1-s;
            
        elseif isnan(nR(i))
            
            if i==na
%                 disp('drive on i')
                Wf(i,j) = W(j)*(nW(j)*(1-hDW*s) +nR(j)*(1-hN*s) +nP(j)*(1-hDP*s))/(nW(j)+nR(j)+nP(j));
            else
%                 disp('N on i')
                Wf(i,j) = W(j)*(1-hN*s);
            end
            
        elseif isnan(nR(j))
            
            if j==na
%                 disp('drive on j')
                Wf(i,j) = W(i)*(nW(i)*(1-hDW*s) +nR(i)*(1-hN*s) +nP(i)*(1-hDP*s))/(nW(i)+nR(i)+nP(i));
            else
%                 disp('N on j')
                Wf(i,j) = W(i)*(1-hN*s);
            end
            
        else
%             disp('i and j functional')
            Wf(i,j) = W(i)*W(j);
        end
        
    end
    
    
end

% Wf





%% Mutation matrix

%For triplex I need to loop over all haplotypes and count the number of W
%R, N and P alleles in each haplotype
% Assume that with partial resistance allele that rate of generation of
% functional is still xi = xiR +xiP
xi = xiR+xiP;

muvec = [1-muSite;xiR*muSite;(1-xi)*muSite;xiP*muSite];
uR = [0;1;0;0];
uN = [0;0;1;0];
uP = [0;0;0;1];

nalleles = zeros(na-1,n-1);

for jj=1:m
    ndims(jj) = n-1;
end
% ndims
M = zeros(na);

for l=1:na-1
%     [i, j, k] = ind2sub(ndims,ind_ud(l));
%     
%     u = [i,j,k];
    
    u = linearind2sub(m,n-1,l);
    
%     % nn is a (n-1)x1 vector of the number of times that allele appears in index u 
%     for kk=1:n-1
%         nalleles(l,kk) = numel(find(u==kk));
%     end
    
    V=1;
    
    for r=1:m
        
        switch u(r)
        
            case 1
                mvec = muvec;        
            
            case 2
                mvec = uR;
        
            case 3
                mvec = uN;
            
            case 4
                mvec = uP;
                        
            
        end
%         mvec
        
        V = kron(mvec,V);
                
            
%          % if allele rr appears nalleles times then multiple mutation
%          % matrix
%          for rr=1:nalleles(l,r) %If nalleles=0 then this should be by-passed
%              V = kron(V,mvec);
%          end
% Order now matters!


        
    end
    
%     V
         
    %Sum over rows of nalleles should always be n-1 and so v should be 3
    %dimensional for triplex and m-dimensional in general - but flattened to 2D
    if m>1
        V = reshape(V,ndims);
    end
    
%     %Do analogue of suming transpose of matrix (V+V^T), of non-diagonal
%     %part only!
%     V = Symmetrise3DTensor(V);
% 
%     v = V(ind_ud);
    v = V(:);
    

    
    M(:,l) = [v;0];
    
    
end

M(na,na) = 1;

Mu = M-eye(na);
% 
% if mu==0 & numel(find(Mu~=0))==0
%     disp('Good')
% else
%     disp('Good')
% end

% Mu

%% Drive heterozygotes

%Set up drive heterozygotes - Drive allele is now last for allele vector

%Now drive only affects gamete production for those heterozygotes with "DDD" 
%that have at least one W in their haplotype, which are 
%WWW/DDD, WWR/DDD, WRR/DDD, WWN/DDD, WRN/DDD, WNN/DDD
%However, approach we'll take is to calculate all combinations using base
%stoichiometry vectors kappa, kR, kN, as below

% Calculate stoichiometry vectors

%This vector kappa is the fraction of gametes expected from W at a single
%site with drive on the other chromosome
%With partial resistance allele P, we now have betaR as fraction of NHEJ
%producing R and betaP fraction producing P, where we can let 

betaW = betaWR+betaWP;
kappaW = [1-epsilonW ; epsilonW*nuW*betaWR;epsilonW*nuW*(1-betaW);epsilonW*nuW*betaWP;epsilonW*(1-nuW)];

%Partial resistance drive "reaction" determines exactly how P alleles when
%paired with D produce different numbers of 
betaP = betaPR+betaPW;
kappaP = [epsilonP*nuP*betaPW ; epsilonP*nuP*betaPR;epsilonP*nuP*(1-betaP);1-epsilonP;epsilonP*(1-nuP)];

%These are the fractions of gametes expected for R and N respectively
kR = [0 1 0 0 0]';
kN = [0 0 1 0 0]';
% kD = [0 0 0 1];

%Here for triplex we need ndims4 = [n,n,n], as well as ndims =
%[n-1,n-1,n-1]
for jj=1:m
    ndims4(jj) = n;
end

S = zeros(na);

%want to loop over all possible haplotypes which could be heterozygote with
%DDD (e.g. XYZ/DDD) where {X,Y,Z} = {W,R,N}, so not including D — since XYD
%is assumed to always be converted to "DDD" and so cannot arise
%So will use linear indexing over na-1 haplotypes not inc D



for l=1:na-1
    
%     VV = zeros(na-1,1);
    
%     [i, j, k] = ind2sub(ndims,ind_ud(l)); %(e.g. l=3-> u=[1,2,2] =WRR
%     
%     u = [i,j,k];
    
    u = linearind2sub(m,n-1,l);
    
%     % nalleles is a (ns-1)x(n-1) vector:nalleles(l,kk) is number of times that allele kk appears in index u/haplotype l (ell not one) 
%     for kk=1:n-1
%         nalleles(l,kk) = numel(find(u==kk));
%     end
    
    %e.g. for l=3-> u=[1,2,2] =WRR, we have nalleles(l,:) = [1,2,0]
    
    V=1;
    for r=1:m
        
        switch u(r)
        
            case 1
                kvec = kappaW;   
            
            case 2
                kvec = kR;
        
            case 3
                kvec = kN;
                
            case 4
                kvec = kappaP;
                
                     
            
        end
                
         
        
%          for rr=1:nalleles(l,r) %If nalleles=0 then this should be by-passed
%              V = kron(V,kvec);
%          end


        V = kron(kvec,V);
         
    end
         

    %Sum over rows of nalleles should always be n-1 and so v should be an 3
    %dimensional for triplex and m-dimensional in general - but flattened to 2D
    
%     V = reshape(V,ndims4);

%Flatten to 1D array — need to check order

    if m>1
        V =reshape(V,numel(V),1);
    end
    
        
%     %Do analogue of suming transpose of matrix (V+V^T), of non-diagonal
%     %part only!
%     V = Symmetrise3DTensor(V); %nxnxn tensor
%     
%     %Drive entry corresponds to the "last" page of the tensor V, and only
%     %the "upper diagonal" (i<=j<=k) 
%     
%     ii = 1:(n)^m;
%     [i, j, k] = ind2sub([n,n,n],ii);
% 
%     ind_ud4 = find(k>=j & j>=i & (i==4 | j==4 | k==4));
%     
%     SD = sum(V(ind_ud4));


%Sum over all indices/haplotypes with at least one D
    SD = sum(V(indDrive));
    
    %Reduce to (n-1)x(n-1)x(n-1) to focus on haplotypes that have only W,R,N
%     V = V(1:n-1,1:n-1,1:n-1);

    V = V(~indDrive);
    
    %ndims defined above for mutation matrix [n-1 n-1 ... n-1] with m
    %elements to represent a m-dimensional array with reduced sub-space of
    %alleles without drive D — so here we reshape V to fit this, and
    %because of linear "column"-wise ordering, we know that taking V without
    %any elements that have at least one D, leaves us with exactly an array
    %of n-1 x n-1 x... n-1
    if m>1
        reshape(V,ndims);
    end
    
    v = V(:);
 
%     %Find "upper diagonal"
%     v = V(ind_ud);
    
    S(:,l) = [v;SD];
    
    
end


S(na,na) = 1;

%Subtracting off identity matrix means the haplotypes corresponding to columns of S which are not
%involved in drive reactions are now zero (e.g. RRN, or RRR, but not WPR,
%since WPR, has a W and a P that can be involved with drive (if in
%heterozygote with DDD)
S = S-eye(na);

%So columns which have non-zero correspond to those haplotypes involved in
%drive reactions
Het = find(sum(abs(S))~=0);

%The complement of this also defines which haplotypes are completely
%resistant e.g. RNR, but not those partially resistant like PRR
%Might want to change this later if effective resistance or population
%elimination can happen even with say PRR
R = setdiff(1:na-1,Het);
test_allele = R;

S = S(:,Het);

Het = [Het;na*ones(size(Het))];

% R
% S
% Het



%% Initial Frequencies
x0 = zeros(na,1);
y0 = zeros(na,1);


%For preexisting=1, resistance alleles are potentially in mutation 
%selection balance — for m>1 the distribution of the frequency of WR, RR is 
%not easy to calculate so run each simulation for time 1/σ to implicitly give
%the frequency distribution of resistance alleles when drive is introduced.
%We can calculate approximately the mean frequency and we use this as the 
%starting frequency of this equilbration phase


if preexisting==1
    
        
    %Calculate approximate mean frequency of resistance alleles in mutation
    %selection balance.
    
    %N.B. That following matrix inversion approach is only valid if all the
    %R containing haplotypes are sufficiently deleterious that they should
    %be at small frequency before drive is introduced, such that non-linear
    %selection terms can be ignored
    %This approach is also only valid as long as the R alleles have no
    %dominance (h=1/2) — although it should be possible to extend this
    %matrix linear approach to include dominance..(?)
    indm = find(Mu(:,1)~=0);
    indm = indm(2:end);
    
    muin = Mu(indm,1);

    f = diag((Wf(1,indm)-1)/2);
    

    MMu = Mu(indm,indm);
    
    x0(indm) = -inv(MMu+f)*muin;
    y0(indm) = -inv(MMu+f)*muin;
    

    %Run simulations for 1/sigma generations *without drive*

    x0(na) = 0;
    x0(1) = 1-sum(x0); 

    y0(na) = 0;
    y0(1) = 1-sum(y0);
    
    Tpre = max(round(1./sigmaR));
    
    [t,x,y,z,M,F,N,wm,wf,resistance,tau] = TargetSiteResistanceSims_while(na,K,psi,x0,y0,Wm,Wf,Mu,S,Het,Rm,Tpre,ratio,Gauss, test_allele, inf, Cswitch);
                                          
    
    x0 = zeros(na,1);
    y0 = zeros(na,1);
    
    x0(indm) = x(indm,end);
    y0(indm) = y(indm,end);
    
    x0(na) = xD;
    x0(1) = 1-sum(x0); 

    y0(na) = yD;
    y0(1) = 1-sum(y0);
    
    
    if ~isempty(find(x0<0, 1))
        disp('Error: frequency vector must not have negative entries')
%         input('')
    end
 
    if ~isempty(find(y0<0, 1))
        disp('Error: frequency vector must not have negative entries')
%         input('')
    end




else
    
    if numel(xinit)>1
        x0 = xinit;
        y0 = yinit;

    else

        xD = xinit;
        yD = yinit;
        
        x0(na) = xD;
        y0(na) = yD;
        
        x0(1) = 1-xD;
        y0(1) = 1-yD;

    end
    
end


% T
[t,x,y,z,M,F,N,wm,wf,resistance,tau] = TargetSiteResistanceSims_while(na,K,psi,x0,y0,Wm,Wf,Mu,S,Het,Rm,T,ratio,Gauss, test_allele, test_threshold, Cswitch);
                                                                     
Nend = N(end);
% tau = t(end);

if plotfig ==1
    
    if vis==0
        f = figure('visible','off')
    else
        figure
    end
    
    col=get(gca,'ColorOrder');
    c = distinguishable_colors(na-7);
  
    col = [col;c];

    subplot(2,1,1)
    
    for j=1:na-1
%         jj=j-1;
        if plotmale==1
            plot(t,x(j,:),'Color',col(j,:), 'LineWidth',1,'LineStyle','-');hold on
        end
        
        plot(t,y(j,:),'Color',col(j,:),'LineWidth',2,'LineStyle','-');hold on
    end
    
    colD = 'k';
    if plotmale==1
        plot(t,x(na,:),'Color',colD, 'LineWidth',1,'LineStyle','-');hold on
    end
    
    plot(t,y(na,:),'Color',colD, 'LineWidth',2,'LineStyle','-');hold on
  
    if plotmale==1
        plot(t(1:end-1),wm(1:end-1),'Color','k','LineWidth',1,'LineStyle','-.');
    end
    plot(t(1:end-1),wf(1:end-1),'Color','k','LineWidth',2,'LineStyle','-.');
    
    
    hold off
   
    xlabel('Time $t$ (generations)')
    ylabel('Frequency')

    pr = 2;
    ht=title(['$m=$',num2str(m,pr)',' gRNAs: $N = ',num2strpow(K,pr),'; s = ',num2str(s),...
        '; h_{DW} = ',num2str(hDW,pr),'; h_{DP} = ',num2str(hDP,pr),'; h_N = ',num2str(hN,pr),...
        '; \sigma_R = ',num2str(sigmaR,pr),'; \sigma_P = ',num2str(sigmaP,pr),...
        '; \nu_W=',num2strpow(nuW,pr),'; \nu_P=',num2strpow(nuP,pr),'; \mu = ',num2strpow(muSite,pr),...
        '; \beta_{WR} = ',num2strpow(betaWR,pr),'; \beta_{WP} = ',num2strpow(betaWP,pr),...
        '; \beta_{PW} = ',num2strpow(betaPW,pr),'; \beta_{PR} = ',num2strpow(betaPR,pr),...
        '; \xi_R = ',num2strpow(xiR,pr),'; \xi_P = ',num2strpow(xiP,pr),'$']);
%     title(['$N = ',num2strpow(K),'$'])
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    set(gca,'FontSize',18)
    set(ht,'FontSize',18)
    
    % Legend
    LS  = KronStrings('WRNP',m);
    LS = vertcat(LS,{'$\mathrm{D}^m$'});

    for k=1:na
        if plotmale==1
            legendstr{2*k-1} = ['Male ',LS{k}];
            legendstr{2*k} = ['Female ',LS{k}];
        else
            legendstr{k} = ['Female ',LS{k}];
        end
            
    end
    
    if plotmale==1
        legendstr{2*na+1} = 'Mean Male Fitness';
        legendstr{2*na+2} = 'Mean Female Fitness';
    else
        legendstr{na+1} = 'Mean Female Fitness';
    end
    
    hleg = legend(legendstr, 'NumColumns',2);
    hleg.FontSize = 12;
    
    
    xlim([0 max(t)])


    subplot(2,1,2);
    plot(t,M,'LineWidth',2);hold on
    plot(t,F,'LineWidth',2);hold on
    plot(t,N,'k','LineWidth',2)
    plot(t,K*ones(size(t)),'k--');hold off
    
    xlabel('Time $t$ (generations)')
    ylabel('Population size')
    
    legend('Male','Female','Total','Intrinsic carrying capacity','Location','southwest')
    
    ylim([0.1 1.2*K])
    xlim([0 max(t)])
    set(gca,'FontSize',22)
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')

%     if vis==0
%         savefig(f,['TargetSiteResistance_Triplex_4Allele_N=',num2strexp(K),'_s=',num2str(s),'_hDW=',num2str(hDW),'_hN=',num2str(hN),'_sigma=',num2str(sigmaR),...
%             '_epsilon=',num2str(epsilonW),'_nu=',num2str(nuW),'_beta=',num2str(betaWR),'_xi=',num2str(xiR),'.fig'])
%     end
    
end

save('Params.mat')

if ~isempty(name)

    mkdir(name)
    cd(name)
    
%     size(t)
%     size(x)
%     size(y)
%     size(z)
%     size(F)
%     size(M)
%     size(N)
%     size(wf)
%     size(wm)

    if m==1
    
        Z = [t',y',x',F',M',N',wf',wm']
        
        results = array2table(Z,...
            'VariableNames',...
            {'time (generations)',...
            'Female W','Female R1','Female R2','Female R3','Female D',...
            'Male W','Male R1','Male R2','Male R3','Male D',...
            'Female Population Size','Male Population Size','Total Population Size',...
            'Mean female fitness', 'Mean male fitness'});
    
        save([name,'.mat'])
        writetable(results,[name,'.csv'])

    end

    if plotfig==1
        
        if vis==0
            savefig(f,[name,'_N=',num2strexp(K),'_s=',num2str(s),'_hDW=',num2str(hDW),'_hN=',num2str(hN),'_sigma=',num2str(sigmaR),...
                '_epsilon=',num2str(epsilonW),'_nu=',num2str(nuW),'_beta=',num2str(betaWR),'_xi=',num2str(xiR),'.fig'])
        else
            savefig([name,'_N=',num2strexp(K),'_s=',num2str(s),'_hDW=',num2str(hDW),'_hN=',num2str(hN),'_sigma=',num2str(sigmaR),...
                '_epsilon=',num2str(epsilonW),'_nu=',num2str(nuW),'_beta=',num2str(betaWR),'_xi=',num2str(xiR),'.fig'])
        end
            
    end
    
    

    cd ..


end



% size(N)
% size(M)
% size(F)
% size(x)
% size(y)
% size(t)
% max(t)
end



function S = SumTensorDiagonal(n,d)

%Number of Alleles n
%Number of dimensions d

S = gamma(n+d)/gamma(n)/factorial(d);



end

function SS = SuperAlleleString(a)

S = OuterProductStrings(a);

SS = [];

for k=1:numel(a)
    for j=1:numel(a)
        SS = horzcat(SS,S(j,k));
    end
end


SS = horzcat(SS,'D');

end


function A = num2strexp(A)

A = num2str(A,'%10.1e');

end




function S = Symmetrise3DTensor(A)


%This produces a symmetrised tensor, where the sum over the 
%"upper diagonal" i<=j<=k, gives the same as the sum over the whole original tensor 
%This is analogous to S = A'+tril(A,-1)' in 2D, i.e. folding back and
%summing the off-diagonal elements leaving diagonal untouched.


%Size of tensor
[N1,N2,N3] = size(A);

if N1~=N2 | N1~=N3 | N2~=N3
    disp('Error: tensor should be "square" (cube)')
else
   nD=N1;
end

%This needs changing to account for different redundancies of permutations

p = perms([1 2 3]);
[N,~] = size(p);

% N
%Extract diagonal
for n=1:nD
    D(n) = A(n,n,n);
end

S = zeros(nD,nD,nD);

%Sum over all permutations
for n=1:N
    S = S + permute(A,p(n,:));
end


%This overcounts for redundant pairs of indices

for i=1:nD
    for j=1:nD
        for k=1:nD
            u = [i,j,k];
            
            %If any two indices are the same, but 3rd not
            
            if i==j & (k~=i & k~=j) | i==k & (j~=i & j~=i) | j==k & (i~=k & i~=k)
                
                S(i,j,k) = S(i,j,k)/2;
            end
        end
    end
end

                




% S = S/6;


for n=1:nD
    S(n,n,n) = D(n);
end

end





function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);

% Copyright 2010-2011 by Timothy E. Holy

  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end

  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end



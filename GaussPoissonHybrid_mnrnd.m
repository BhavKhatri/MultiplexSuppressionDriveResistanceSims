function n = GaussPoissonHybrid_mnrnd(N,x)

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

%Hybrid Gaussian-Poisson approx to multinomial random number generator
%N should be an integer
%x should be a vector of frequencies/probabilities and whose elements sum
%to one
%n is an output vector of same length as x and whose elements sum to N


if ~isempty(find(x<0, 1))
    disp('Error: frequency vector must not have negative entries')
end

if abs(sum(x)-1)>3*eps
    disp('Error: elements of frequency vector must sum to 1')
end


%Order frequencies from small to large and record ordering
%This is needed as mvrnd seems to have problems when there is a large
%dynamic range along diagonal of covariance matrix
[x, indsort] = sort(x);

%Find reverse indexing
unsorted = 1:length(x);
reverse_ind(indsort) = unsorted;


if N<1e3
    %The approximation requires that the threshold used nth<<N — nth will
    %typically be about 10 so if N<1e3 seems reasonable to just use exact
    %multinomial random numbers — for N<1e3 multinomial in matlab runs fast
    n = mnrnd(N,x);
    n = n';
else

    %Threshold to use Poisson appprox
    nthr = 10;


    indnot0 = x~=0;



    n = zeros(size(x));

    %find all entries which have expected number in next generation <=gamma
    indsmall = N*x <= nthr;
%     N*x
    
    if isempty(find(indsmall==1, 1)) %There are no "small" entries 
        
        nsum = 0;
        
    else

        %Instead treating these "multinomially" as they are small compared to
        %N:
        %draw from Poisson distribution instead
        n(indsmall & indnot0) = poissrnd(N*x(indsmall & indnot0));

        %Sum the number of "small" alleles
        nsum = sum(n((indsmall & indnot0)));
        

        %index ~indsmall are those alleles with expected number >nthr & by def
        %~=0
    
    end

    NN = N-nsum;
    
    
    %Convert frequencies of "large" alleles to frequencies within subset of
    %large alleles and also take only on simplex

    indlarge  = find(N*x > nthr);

    if numel(indlarge)>1 %more than one large allele & there has to be at least 1!

        %xx should be on the simplex of the subset of alleles which
        %are "large"
        xx = x(indlarge(1:end-1))/sum(x(indlarge));
        %So xx does not sum to one


        %Draw random number from multivariate normal distribution with
        %variance matrix of allele frequencies given by B
% 
%         B = AlleleVarianceMatrix(xx)/NN
        B = NN*AlleleVarianceMatrix(xx);


%         nn = NN*mvnrnd(xx,B)
        nn = mvnrnd(NN*xx,B);



    
        while sum(nn)>NN || ~isempty(find(nn<0, 1)) %So rarely negative entries can be produced
%             nn = NN*mvnrnd(xx,B);
            nn = mvnrnd(NN*xx,B);
        end
        
        
        %The Gaussian approx will produce non-integer nn for the large
        %subset=> round them
        nn = round(nn); 
        n(indlarge(1:end-1)) = nn;         
          
        %The total sum of the large subset should equal NN
        nlast = NN - sum(nn);
        n(indlarge(end)) = nlast;
        
        
    else
        %There is only a single large allele and it's number must be NN
        n(indlarge)=NN;
        

    end


   
end

%Revert back to original order of indexing
n = n(reverse_ind);

end


function B = AlleleVarianceMatrix(x)


B = zeros(numel(x));


for j = 1:numel(x)
    for k=1:numel(x)
        
        if j==k
            B(j,k) = x(j)*(1-x(j));
        elseif k>j %Only fill in upper triangular part
            B(j,k) = -x(j)*x(k);
        end
        
    end
end
        
        
B = triu(B,1)'+B;

end
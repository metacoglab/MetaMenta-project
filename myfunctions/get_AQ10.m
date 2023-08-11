function [AQ10, AQ10_AD, AQ10_AS, AQ10_C, AQ10_S, AQ10_I] = get_AQ10(ratings)

%ratings = participants' ratings on Autism Quotient (1-10 from 1
%'Definitely Agree' to 4 'Definitely disagree'
%nsubjects = number of subjects (also the number of rows of the ratings)

%returns AQ10 rating 
%elisavanderplasATgmail.com


AQ_AD = 0;
AQ_AS = 0; 
AQ_C = 0; 
AQ_S = 0; 
AQ_I = 0;
    
for i = [4,10,7,3,5,6] %reverse keyed items, note these are the order of the items as presented in Gorilla, not the order in the original questionairre 
    ratings(i) = 5-ratings(i);
end

    AQ10 = sum(ratings>2); %all items
    
    AQ_AD = [AQ_AD, ratings(1), ratings(4)]; 
    AQ_AS = [AQ_AS, ratings(7), ratings(10)];
    AQ_C = [AQ_C, ratings(3), ratings(5)];
    AQ_I = [AQ_I, ratings(2), ratings(8)];
    AQ_S = [AQ_S, ratings(6), ratings(9)]; 
        
    AQ10 = sum(ratings>2); %all items
    
    AQ10_AD = sum(AQ_AD>2); 
    AQ10_AS = sum(AQ_AS>2); 
    AQ10_C = sum(AQ_C>2); 
    AQ10_S = sum(AQ_S>2); 
    AQ10_I = sum(AQ_I>2); 
end

function [PC]=Hermite_polynomial_generation(dim_max, order_max)

%% set required dimension and order

% clc;clear;
% dim_max=1;
% order_max=8;



%% iterative algorithm to generate all the combinations
% the algorithm here follows a classical programming exercise "Letter
% Combinations of a Phone Number" in LeetCode. Iterative approach is
% employed here, and recurve approach which can also solve the problem is not explored here.


% for dim=1:dim_max
%    if dim==1
%         results=zeros(order_max+1,dim_max);
%         for i=0:order_max
%             results(i+1,1)=i;
%         end
%         continue;
%    end
% 
%    old_mat_size=size(results,1);
%    new_mat_size=(order_max+1)*old_mat_size;
%    temp_results=zeros(new_mat_size, dim_max);
%    count=1;
%    for order=0:order_max
%        for i=1:old_mat_size
%            temp_PC=results(i,:);
%            temp_PC(dim)=order;
%            if sum(temp_PC)<=order_max
%                temp_results(count,:)=temp_PC;
%                count=count+1;
%            end
%            
%        end
%   
%    end
%    temp_results=temp_results(1:count-1,:);
%    results=temp_results;
%    
%     dim
% end

%% iterative algorithm to generate all the combinations

for dim=1:dim_max
   if dim==1
        results=zeros(order_max+1,dim_max);
        for i=0:order_max
            results(i+1,1)=i;
        end
        continue;
   end

   old_mat_size=size(results,1);
   new_mat_size=(order_max+1)*old_mat_size;
   temp_results=sparse(new_mat_size, dim_max);
   count=1;
   for order=0:order_max
       results_with_order=results;
       results_with_order(:,dim)=order;
       temp_results(order*old_mat_size+1:(order+1)*old_mat_size, :)=results_with_order;
   end
   temp_sum=sum(temp_results,2);
   index=(temp_sum<=order_max);
   temp_results=temp_results(index,:);
   results=temp_results;
   
    dim
end


%%
sum_results=sum(results,2);
[sum_results_sorted,I] = sort(sum_results, 'ascend');

final_results_sorted=results(I,:);
PC=final_results_sorted;




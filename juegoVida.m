clc;
clear;
close all;


N = input("Tamano N: ");
num = input("Numero de generaciones: "); % pasos

A = rand(N);
p= 0.3; % probabilidad de vida de la primera generacion

% estado inicial
for (i=1:N)
    for (j=1:N)
        if A(i,j)<p
            A(i,j)=true;
        else
            A(i,j)=false;
        end
    end
end

% contorno
B = false(N+2,N+2,num+1);
B(2:N+1,2:N+1,1)=A;

% estado inicial
spy(B(:,:,1));

for i=2:num+1
    for j=2:N+1
        for k=2:N+1
            M = B(j-1:j+1,k-1:k+1,i-1);
            % contamos los vecinos con vida
            sum_ones = sum(M(:));
            if (B(j,k,i-1) == true)
                if(sum_ones==3 || sum_ones==4)
                    B(j,k,i)=true; % se mantiene con vida
                end
                %en otro caso la celula se muere
            else
                if (sum_ones==3)
                    B(j,k,i)=true; % revive
                end
            end
        end
    end
    pause(0.2);
    spy(B(:,:,i)); % plot de la generacion i
    drawnow;

end

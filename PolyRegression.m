% This script applies Newton's Method to get optimized solution for poly-regress

% You can chagne "d" to test different polynomial powers (default is cubic polynomial)
d = 4;

% You can chagne "Xs" and "Ys" to test different datasets
Xs = [-3.00 -2.67 -2.33 -2.00 -1.67 -1.33 -1.00 -0.67 -0.33 0.00 0.33 0.67 1.00 1.33 1.67 2.00];
Ys = [4.54 1.94 -0.72 -1.22 -2.47 -2.61 -2.63 -2.38 -2.39 -1.61 -1.80 -1.35 -1.17 -1.44 -1.66 -2.59];
N = length(Xs);

figure();
hold on;
	scatter(Xs,Ys,'filled');
	ps = [1.5,2.0,5.0]; colors = ['b-','g-','r-'];
	for i = 1 : 3
		qInitial = randn(d,1); 
		[q,~,fe] = newton(qInitial,Xs,Ys,ps(i));
		DATA_i = linspace(min(Xs),max(Xs),4*N);
		plot(DATA_i,polyval(q(end:-1:1),DATA_i), colors(i));
		fprintf(colors(i)+' p=%0.2f, y=%0.2fx^3+(%0.2f)x^2+(%0.2f)x+(%0.2f), err=%0.2f', ps(i), q(end:-1:1), fe);
    	end
hold off;

% Gets the p-error of the simulating polynomial
function result = Ep(q,x,y,p)
    result = sum(abs(polyval(q(end:-1:1),x)-y).^p);
end

% Gets the first derivative of Ep
function result = dEp(q,x,y,p)
    r = polyval(q(end:-1:1),x)-y;
    result = zeros(length(q),1);
    for j=0:length(q)-1
        result(1+j) = sum(p*(abs(r).^(p-1)).*(x.^j).*sign(r));
    end
end

% Gets the second derivates of Ep
function result = ddEp(q,x,y,p)
    r = polyval(q(end:-1:1),x)-y;
    result = zeros(length(q),length(q));
    for j=0:length(q)-1
        for k=0:length(q)-1
            result(1+j,1+k) = sum(p*(p-1)*(abs(r).^(p-2)).*(x.^(j+k)));
        end
    end
end

% Gets optimal polynomial minimizing Ep
function [root,newtonErr,err] = newton(qInitial,x,y,p)
    maxItr = 15; itr = 0; q = qInitial;
    while (itr < maxItr)
        dq = ddEp(q,x,y,p)\dEp(q,x,y,p);
        q = q - dq;
        itr = itr+1;
    end
    root = q; newtonErr = norm(dEp(root,x,y,p)); err = Ep(root,x,y,p);
end

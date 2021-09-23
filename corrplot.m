% plot two variables and display their correlation coefficient

% Arun Sripati
% January 13 2009

function [c,p,hobj] = corrplot(x,y,titlestr,unityflag,linespec,ctype)
if(~exist('titlestr')|| isempty(titlestr)),titlestr = ''; end 
if(~exist('linespec') || isempty(linespec)),linespec = '.'; end 
if(~exist('unityflag')|| isempty(unityflag)),unityflag = 0; end 
if(~exist('ctype') || isempty(ctype)),ctype = 'pearson'; end 

x = x(:); y = y(:); 

hobj=plot(x,y,linespec); c = get(hobj,'Color'); hold on; 
q = find(isfinite(x.*y)); a = polyfit(x(q),y(q),1); yp = polyval(a,x); plot(x,yp,'Color',c); 
axis square; 
xlabel(inputname(1)); ylabel(inputname(2)); 
if(unityflag)
    m = [min(min(x),min(y)) max(max(x),max(y))]; axis([m m]); 
    unityslope(1,'k');
end
% [c,p] = nancorrcoef(x,y);
[c,p] = corr(x,y,'rows','complete','type',ctype);
n = length(find(isfinite(x)&isfinite(y))); 
if(isempty(titlestr)),
    title(sprintf('n = %d, r = %2.2f, p = %2.2g',n,c,p));
else
    title(sprintf('%s (n = %d, r = %2.2f, p = %2.2g)',titlestr,n,c,p));
end

return

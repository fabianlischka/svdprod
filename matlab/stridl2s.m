function TS = stridl2s( TL )
% STRIDL2S turns a symmetric tridiagonal matrix TL in usual ("long") format
% into "short" format, with TS(:,1) being the diagonal, and TS(:,2) the
% subdiagonal (and TS(end,2) being zero)

TS = [ diag(TL) [ diag(TL,1); 0 ] ];
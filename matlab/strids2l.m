function TL = strids2l( TS )
% STRIDS2L turns a symmetric tridiagonal matrix TS in "short" format
% into normal "long" format, with TS(:,1) being the diagonal, and TS(:,2) the
% subdiagonal (and TS(end,2) being zero)

TL = diag(TS(:,1)) + diag(TS(1:end-1,2), 1) + diag(TS(1:end-1,2), -1);
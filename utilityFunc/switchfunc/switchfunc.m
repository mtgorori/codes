function [ out ] = switchfunc( varargin)
%SWITCHFUNC  条件分岐を簡潔に書くための関数
%   SWITCHFUNC( cond_1, value_1,...
%                        cond_2, value_2,...
%                        true, default_value)
%cond_1, cond_2, の中で最初にtrueの値をとるものを
%cond_kとすると，value_kを値として返す．
%全ての条件が満たされない場合，エラーを返す．
for k=1:2:nargin-1
    if varargin{k}
        out = varargin{k+1};
        return
    end
end
error('switchfunc: nomatch');
end


function [ out ] = switchfunc( varargin)
%SWITCHFUNC  ����������Ȍ��ɏ������߂̊֐�
%   SWITCHFUNC( cond_1, value_1,...
%                        cond_2, value_2,...
%                        true, default_value)
%cond_1, cond_2, �̒��ōŏ���true�̒l���Ƃ���̂�
%cond_k�Ƃ���ƁCvalue_k��l�Ƃ��ĕԂ��D
%�S�Ă̏�������������Ȃ��ꍇ�C�G���[��Ԃ��D
for k=1:2:nargin-1
    if varargin{k}
        out = varargin{k+1};
        return
    end
end
error('switchfunc: nomatch');
end


function y = logger(formatSpec, varargin)
% ���O�p�̊֐�: �W���o�͂���у��O�t�@�C���ւ̏o�͂��s��
% ��{�I�ɂ͕W���o�͂�fprintf�Ɠ����悤�Ɏg���΂悢
% ���O�ɂ͌Ăяo�����̊֐����ƌĂяo�����s��������
% ���̊֐��͌Ăяo�����Ƀt�@�C�����J���Ă���

% formatSpec: ����������
% varargin: �����ɓ��ꍞ�ޕϐ��i�����j

st = dbstack(); % �֐��X�^�b�N�Ɋւ���\����
callfname = st(2).name; % �Ăяo�����̊֐��̖��O
callline = num2str(st(2).line); % �Ăяo�����֐��ł̌Ăяo���s

% �֐��̖��O�ƌĂяo���s�͕�����ŗ^������
% ���������O��ŕ��ׂ��Ƃ��C���ꂢ�Ɍ�����悤�ɊԊu����������
len1 = strlength(callfname); 
len2 = strlength(callline);
maxlen1 = 17;
maxlen2 = 4;
bl1 = blanks(abs(maxlen1 - len1)); % �֐����̌�̋󔒂̐�
bl2 = blanks(abs(maxlen2-len2)); % �s�����̌�̋󔒂̐�

% �󔒂ɂ��e������̊Ԋu����
formatSpec2 = ['$ [' callfname '] line:' callline bl1 bl2 '| ' formatSpec]; 

%-----------------
% �W���o�͂ɏo��
%-----------------

fprintf(formatSpec2, varargin{1:nargin-1}); 

%------------------
% �t�@�C���ւ̏o��
%------------------

fid = fopen('./logger.txt', 'a'); % �t�@�C�����J��

fprintf(fid, formatSpec2, varargin{1:nargin-1}); % �t�@�C���ɏ�������
fclose(fid); % �t�@�C�������

y = 0;
end
function speechproc()
    % ���峣��
    FL = 80;                % ֡��
    WL = 240;               % ����
    P = 10;                 % Ԥ��ϵ������
    Fs = 8000;              % ������
    s = readspeech('voice.pcm', 100000);            % ��������s
    L = length(s);          % ������������
    FN = floor(L / FL) - 2;     % ����֡��
    % Ԥ����ؽ��˲���
    exc = zeros(L, 1);       % �����źţ�Ԥ����
    zi_pre = zeros(P, 1);    % Ԥ���˲�����״̬
    s_rec = zeros(L, 1);     % �ؽ�����
    zi_rec = zeros(P, 1);    % �ؽ��˲�����״̬
    % �ϳ��˲���
    exc_syn = zeros(L, 1);   % �ϳɵļ����źţ����崮��
    zi_syn = zeros(P, 1);    % �ϳ��˲�����״̬
    s_syn = zeros(L, 1);     % �ϳ�����
    % ����������˲���
    exc_syn_t = zeros(L, 1);   % �ϳɵļ����źţ����崮��
    zi_syn_t = zeros(P, 1);    % �ϳ��˲�����״̬
    s_syn_t = zeros(L, 1);     % �ϳ�����

    exc_syn_tL1 = zeros(L, 1);
    zi_syn_tL1 = zeros(P, 1);
    s_syn_tL1 = zeros(L, 1);
    exc_syn_tL2 = zeros(L, 1);
    zi_syn_tL2 = zeros(P, 1);
    s_syn_tL2 = zeros(L, 1);
    exc_syn_tS1 = zeros(L, 1);
    zi_syn_tS1 = zeros(P, 1);
    s_syn_tS1 = zeros(L, 1);
    exc_syn_tS2 = zeros(L, 1);
    zi_syn_tS2 = zeros(P, 1);
    s_syn_tS2 = zeros(L, 1);

    % ���ٲ�����˲����������ٶȼ���һ����
    exc_syn_v = zeros(2 * L, 1);   % �ϳɵļ����źţ����崮��
    zi_syn_v = zeros(P, 1);        % �ϳ��˲�����״̬
    s_syn_v = zeros(2 * L, 1);     % �ϳ�����

    hw = hamming(WL);       % ������
    k = 2 * FL + 1;
    k_t = 2 * FL + 1;
    FL_v = 2 * FL;
    k_v = 2 * FL + 1;
    % ���δ���ÿ֡����
    for n = 3: FN
        % ����Ԥ��ϵ��������Ҫ���գ�
        s_w = s(n * FL - WL + 1: n * FL).*hw;    % ��������Ȩ�������
        [A, E] = lpc(s_w, P);            % ������Ԥ�ⷨ����P��Ԥ��ϵ��
        newA = rotAngle(A, 2 * pi * 150 / Fs);
        newAL1 = rotAngle(1.1 * A, 2 * pi * 150 / Fs);
        newAL2 = rotAngle(1.2 * A, 2 * pi * 150 / Fs);
        newAS1 = rotAngle(0.9 * A, 2 * pi * 150 / Fs);
        newAS2 = rotAngle(0.8 * A, 2 * pi * 150 / Fs);
                                        % A��Ԥ��ϵ����E�ᱻ��������ϳɼ���������
        if n == 27
        % �۲�Ԥ��ϵͳ���㼫��ͼ
            disp(A)
            zplane(-A, 1)
        end

        % ��֡�����������Ҫ����������
        s_f = s((n - 1) * FL + 1: n * FL);

        % (4) �ڴ�λ��д������filter����s_f���㼤����ע�Ᵽ���˲���״̬
        [exc((n - 1) * FL + 1: n * FL), zi_pre] = filter(-A, 1, s_f, zi_pre);

        % (5) �ڴ�λ��д������filter������exc�ؽ�������ע�Ᵽ���˲���״̬
        [s_rec((n - 1) * FL + 1: n * FL), zi_rec] = ...
            filter(1, -A, exc((n - 1) * FL + 1: n * FL), zi_rec);
        s_Pitch = exc(n * FL - 222: n * FL);
        PT = findpitch(s_Pitch);    % �����������PT����Ҫ�����գ�
        G = sqrt(E * PT);           % ����ϳɼ���������G����Ҫ�����գ�

        % (10) �ڴ�λ��д�������ɺϳɼ��������ü�����filter���������ϳ�����
        while k <= n * FL
            exc_syn(k) = G;
            k = k + PT;
        end
        [s_syn((n - 1) * FL + 1: n * FL), zi_syn] = ...
            filter(1, -A, exc_syn((n - 1) * FL + 1: n * FL), zi_syn);

        % (11) ���ı�������ں�Ԥ��ϵ�������ϳɼ����ĳ�������һ��������Ϊfilter
        % ������õ��µĺϳ���������һ���ǲ����ٶȱ����ˣ�������û�б䡣
        while k_v <= n * FL_v
            exc_syn_v(k_v) = G;
            k_v = k_v + PT;
        end
        [s_syn_v((n - 1) * FL_v + 1: n * FL_v), zi_syn_v] = ...
            filter(1, -A, exc_syn_v((n - 1) * FL_v + 1: n * FL_v), zi_syn_v);
        
        % (13) ���������ڼ�Сһ�룬�������Ƶ������150Hz�����ºϳ�������������ɶ���ܡ�
        while k_t <= n * FL
            exc_syn_t(k_t) = G;
            exc_syn_tL1(k_t) = G;
            exc_syn_tL2(k_t) = G;
            exc_syn_tS1(k_t) = G;
            exc_syn_tS2(k_t) = G;
            k_t = k_t + round(PT / 2);
        end
        [s_syn_t((n - 1) * FL + 1: n * FL), zi_syn_t] = ...
            filter(1, -newA, exc_syn_t((n - 1) * FL + 1: n * FL), zi_syn_t);

        [s_syn_tL1((n - 1) * FL + 1: n * FL), zi_syn_tL1] = ...
            filter(1, -newAL1, exc_syn_tL1((n - 1) * FL + 1: n * FL), zi_syn_tL1);

        [s_syn_tL2((n - 1) * FL + 1: n * FL), zi_syn_tL2] = ...
            filter(1, -newAL2, exc_syn_tL2((n - 1) * FL + 1: n * FL), zi_syn_tL2);

        [s_syn_tS1((n - 1) * FL + 1: n * FL), zi_syn_tS1] = ...
            filter(1, -newAS1, exc_syn_tS1((n - 1) * FL + 1: n * FL), zi_syn_tS1);

        [s_syn_tS2((n - 1) * FL + 1: n * FL), zi_syn_tS2] = ...
            filter(1, -newAS2, exc_syn_tS2((n - 1) * FL + 1: n * FL), zi_syn_tS2);
    end

    % (6) �ڴ�λ��д������һ�� s ��exc �� s_rec �к����𣬽�����������
    figure
    subplot(3, 1, 1);
    plot((0: L - 1)./Fs, s);
    xlabel("Length");
    ylabel("$s(n)$", Interpreter="latex");
    subplot(3, 1, 2);
    plot((0: L - 1)./Fs, exc);
    xlabel("Length");
    ylabel("$e(n)$", Interpreter="latex");
    subplot(3, 1, 3);
    plot((0: L - 1)./Fs, s_rec);
    xlabel("Length");
    ylabel("$\hat{s}(n)$", Interpreter="latex");

    x = 2200: 3800;
    figure
    subplot(3, 1, 1);
    plot(x./Fs, s(x));
    xlabel("Length");
    ylabel("$s(n)$", Interpreter="latex");
    subplot(3, 1, 2);
    plot(x./Fs, exc(x));
    xlabel("Length");
    ylabel("$e(n)$", Interpreter="latex");
    subplot(3, 1, 3);
    plot(x./Fs, s_rec(x));
    xlabel("Length");
    ylabel("$\hat{s}(n)$", Interpreter="latex");

    s_fft = fft(s);
    exc_fft = fft(exc);
    s_rec_fft = fft(s_rec);
    % �˴�ͨ��һ���ְ�ѧ���õ���ȷ��Ƶ��ͼ����
    freq = Fs * (0: (L / 2)) / L;
    freqv = Fs * (0: L) / (2 * L);
    figure;
    subplot(3, 1, 1);
    plot(freq, abs(s_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s(n)$", Interpreter="latex");
    
    subplot(3, 1, 2);
    plot(freq, abs(exc_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$e(n)$", Interpreter="latex");

    subplot(3, 1, 3);
    plot(freq, abs(s_rec_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$\hat{s}(n)$", Interpreter="latex");

    figure
    subplot(2, 1, 1);
    plot((0: L - 1)./Fs, s);
    xlabel("Length");
    ylabel("$s(n)$", Interpreter="latex");
    subplot(2, 1, 2);
    plot((0: L - 1)./Fs, s_syn);
    xlabel("Length");
    ylabel("$s_{syn}(n)$", Interpreter="latex");

    s_syn_fft = fft(s_syn);
    figure;
    subplot(2, 1, 1);
    plot(freq, abs(s_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s(n)$", Interpreter="latex");
    
    subplot(2, 1, 2);
    plot(freq, abs(s_syn_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{syn}(n)$", Interpreter="latex");

    s_syn_t_fft = fft(s_syn_t);
    s_syn_tL1_fft = fft(s_syn_tL1);
    s_syn_tL2_fft = fft(s_syn_tL2);
    s_syn_tS1_fft = fft(s_syn_tS1);
    s_syn_tS2_fft = fft(s_syn_tS2);
    figure
    subplot(5, 1, 1);
    plot((0: L - 1)./Fs, s_syn_t);
    xlabel("Length");
    ylabel("$s_{synt}$", Interpreter="latex");
    subplot(5, 1, 2);
    plot((0: L - 1)./Fs, s_syn_tL1);
    xlabel("Length");
    ylabel("$s_{syntL1}(n)$", Interpreter="latex");
    subplot(5, 1, 3);
    plot((0: L - 1)./Fs, s_syn_tL2);
    xlabel("Length");
    ylabel("$s_{syntL2}(n)$", Interpreter="latex");
    subplot(5, 1, 4);
    plot((0: L - 1)./Fs, s_syn_tS1);
    xlabel("Length");
    ylabel("$s_{syntS1}(n)$", Interpreter="latex");
    subplot(5, 1, 5);
    plot((0: L - 1)./Fs, s_syn_tS2);
    xlabel("Length");
    ylabel("$s_{syntS2}(n)$", Interpreter="latex");

    figure
    subplot(5, 1, 1);
    plot(freq, abs(s_syn_t_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{synt}$", Interpreter="latex");
    subplot(5, 1, 2);
    plot(freq, abs(s_syn_tL1_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{syntL1}(n)$", Interpreter="latex");
    subplot(5, 1, 3);
    plot(freq, abs(s_syn_tL2_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{syntL2}(n)$", Interpreter="latex");
    subplot(5, 1, 4);
    plot(freq, abs(s_syn_tS1_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{syntS1}(n)$", Interpreter="latex");
    subplot(5, 1, 5);
    plot(freq, abs(s_syn_tS2_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{syntS2}(n)$", Interpreter="latex");

    s_syn_v_fft = fft(s_syn_v);
    figure
    subplot(3, 1, 1);
    plot((0: L - 1)./Fs, s);
    xlabel("Length");
    ylabel("$s(n)$", Interpreter="latex");
    subplot(3, 1, 2);
    plot((0: 2 * L - 1)./Fs, s_syn_v);
    xlabel("Length");
    ylabel("$s_{synv}(n)$", Interpreter="latex");
    subplot(3, 1, 3);
    plot((0: L - 1)./Fs, s_syn_t);
    xlabel("Length");
    ylabel("$s_{synt}(n)$", Interpreter="latex");

    figure
    subplot(3, 1, 1);
    plot(freq, abs(s_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s(n)$", Interpreter="latex");
    subplot(3, 1, 2);
    plot(freqv, abs(s_syn_v_fft(1: L + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{synv}(n)$", Interpreter="latex");
    subplot(3, 1, 3);
    plot(freq, abs(s_syn_t_fft(1: L / 2 + 1)));
    xlabel("Frequence/Hz");
    ylabel("$s_{synt}(n)$", Interpreter="latex");

    sound(s);
    sound(exc);
    sound(s_rec);
    sound(s_syn);
    sound(s_syn_v);
    sound(s_syn_t);
    sound(s_syn_tL1);
    sound(s_syn_tL2);
    sound(s_syn_tS1);
    sound(s_syn_tS2);

    % ���������ļ�
    writespeech('exc.pcm', exc);
    writespeech('rec.pcm', s_rec);
    writespeech('exc_syn.pcm', exc_syn);
    writespeech('syn.pcm', s_syn);
    writespeech('exc_syn_t.pcm', exc_syn_t);
    writespeech('syn_t.pcm', s_syn_t);
    writespeech('exc_syn_v.pcm', exc_syn_v);
    writespeech('syn_v.pcm', s_syn_v);
return

% ��PCM�ļ��ж�������
function s = readspeech(filename, L)
    fid = fopen(filename, 'r');
    s = fread(fid, L, 'int16');
    fclose(fid);
return

% д������PCM�ļ���
function writespeech(filename, s)
    fid = fopen(filename, 'w');
    fwrite(fid, s, 'int16');
    fclose(fid);
return

% ����һ�������Ļ������ڣ���Ҫ������
function PT = findpitch(s)
    [B, A] = butter(5, 700 / 4000);
    s = filter(B, A, s);
    R = zeros(143, 1);
    for k=1: 143
        R(k) = s(144: 223)'*s(144 - k: 223 - k);
    end
    [R1, T1] = max(R(80: 143));
    T1 = T1 + 79;
    R1 = R1 / (norm(s(144 - T1: 223 - T1)) + 1);
    [R2, T2] = max(R(40: 79));
    T2 = T2 + 39;
    R2 = R2 / (norm(s(144 - T2: 223 - T2)) + 1);
    [R3, T3] = max(R(20: 39));
    T3 = T3 + 19;
    R3 = R3 / (norm(s(144 - T3: 223 - T3)) + 1);
    Top = T1;
    Rop = R1;
    if R2 >= 0.85 * Rop
        Rop = R2;
        Top = T2;
    end
    if R3 > 0.85 * Rop
        Rop = R3;
        Top = T3;
    end
    PT = Top;
return

function ang = rotAngle(v, deltaAngle)
    p = roots(v);
    ang = poly(arrayfun(@(p) p * ...
        exp(deltaAngle * sign(imag(p)) * 1i), p)) * v(1);
return
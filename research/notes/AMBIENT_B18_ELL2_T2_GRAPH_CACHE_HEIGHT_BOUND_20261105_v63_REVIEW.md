# cache-only値job前のLCM・絶対分子上界の独立検証

各既知source classについてD=全分母LCM、U=max|D*c|を独立再計算する。積はDとUの積に最大項数を掛け、和は分母LCMへ正に拡大して絶対分子を加える。全保存envelopeを再構成し、source-derived項数も確認する。原点零・感度零で証明された項だけを省略する。

primitive nu=t+s_trace²/4とraw33積、bar-rho dot、全46 A23、normalized graphまでの保存上界は最大538bitである。factored workerが途中で作るOmegaとS1 pivot補正にも別の正の上界を追加し、2048bit以内を確認する。個別native ell firstにはL0Se+LeS0を含め、effective firstだけでは済ませない。

signed Phi、A23、二不変量、graphの和は一つも評価しない。これは親が明示した今回の2048bit guardの事前証明であり、以前の384bit capを自動変更するものではない。

class info:
    def __init__(self,
                 T1_generate=[1500, 1000],  # vassel, muscle
                 T2=[50, 50],
                 TR=2.8,
                 TI_5=120,
                 TI_3=320,
                 HR=60,
                 # TFE = 49,
                 fa=30,
                 z0=32,
                 bandwidth=10,
                 index_list=[0, 5, 1, 6, 2, 7, 3, 4],
                 thickness=5,
                 b0=1.5,  # Tesla
                 prep_num=1,
                 # tau_x = 1.0, # 原本应该是 receiver bandwidth， 但因为序列设计
                 # receiver_bandwidth = 83.3, # khz
                 gamma=4258,  # Hz / G
                 flow_speed=10,  # cm/s
                 delta=0.1,  # spatial resolution delta
                 fov=12.8,  # cm
                 # tau_y = 0.29

                 ) -> None:
        '''
        在这里，整体的时间如下(按顺序)：
        total_time = (180y + TI_5[0] + TR * rep_time + TI_5[1] - ...)
        其中 TE 为每个 TR 中 10y时刻 到 readout 梯度的中点
        '''

        self.T1 = T1_generate  # float64
        self.T2 = T2
        self.TI_5 = TI_5
        self.TI_3 = TI_3
        self.TE = TR / 2
        # self.TFE = TFE
        self.TR = TR
        self.fov = fov
        self.gamma = gamma
        self.index_list = index_list
        self.z0 = z0
        self.thickness = thickness
        # self.Gx = Gx

        self.b0 = b0
        self.fa = fa
        self.HR = HR
        self.delta = delta

        self.flow_speed = flow_speed
        # TODO: 根据 flow_speed 算出 each_time：相隔多少毫秒，让 flow_speed 表示为每移动1格需要的时间 (ms)
        self.each_time_flow = self.delta / (self.flow_speed * 1e-3)

        # self.w0 = self.gamma * 1e-4 * self.b0 # MHz
        self.w0 = 0
        self.N_pe = int(self.fov / self.delta)
        self.N_read = int(self.fov / self.delta)
        self.prep_num = prep_num

        self.delta_k = 1 / (self.fov)

        self.length = self.N_pe
        self.bandwidth = bandwidth

        # self.BW = receiver_bandwidth # khz
        # self.delta_t = self.tau_x / self.N_read
        # self.delta_t = 1 / (2 * self.BW) # sampling period, ms
        # self.tau_x = self.N_read * self.delta_t
        # self.Gx = 2 * self.BW * 1e3 / (self.gamma * self.fov)

        # self.tau_y = tau_y # ms
        self.tau_y = 1/4 * (self.TR)
        self.tau_x = 2 * self.tau_y
        self.delta_t = self.tau_x / self.N_read

        self.Gx = 1e3 / (self.gamma * self.delta_t * self.fov)
        self.delta_ky = 1 / (self.fov)  # cm-1
        # self.k_max = 1 / (2 * self.delta)
        self.ky_max = 0.5 * self.N_pe * self.delta_ky
        # self.ky_list = [((self.N_pe - 1) / 2 - m) * self.delta_ky for m in np.arange(0, self.N_pe, 1)]
        self.Gyp = self.ky_max * 1e3 / (self.gamma * self.tau_y)
        self.Gyi = self.delta_ky * 1e3 / (self.gamma * self.tau_y)
        # self.m0 = torch.Tensor([0,0,1]).to(device).T

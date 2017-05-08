function Kalman = fill_Kalman(state,Kalman)
            Kalman.kalman = Create_New_Kalman(state.kalman_state,8 ^2*eye(10),2*eye(10));
            Kalman.prev_state = state.kalman_state';
            Kalman.state = state.kalman_state';
            Kalman.BW = state.BW;
            Kalman.Contour = state.Contour;
            Kalman.size = state.size;
            Kalman.weightedSize = state.weightedSize;
            

end

#include<iostream>
#include<cmath>
#include<memory>
#include<vector>

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>

template<typename T>
T pow2(const T num)
{
    return num * num;
}

double normalizeAngle(const double angle)
{
    double res = std::fmod(angle + M_PI, 2*M_PI);
    res = res < 0 ? res + 2*M_PI : res;
    return res - M_PI;
}

struct point
{
    double x, y, theta, v, w;
    point():x(0),y(0),theta(0), v(0), w(0){}
    point(double x_, double y_, double theta_, double v_, double w_):x(x_), y(y_), theta(theta_), v(v_), w(w_){}
    bool operator == (point point2)
    {
        return fabs(x - point2.x) < 1e-4 && fabs(y - point2.y) < 1e-4 && fabs(normalizeAngle(theta - point2.theta)) < M_1_PI;
    }

    bool operator != (point point2)
    {
        return fabs(x - point2.x) >=1e-2 || fabs(y - point2.y) >= 1e-2;
    }
};

struct control
{
    double v, w;
    control():v(0), w(0){}
};


class lqr
{
public:
    void constructMatrix()
    {
        Eigen::Matrix<double, 3, 3> A;
        A<<0, 0, -target_state_.v*std::sin(target_state_.theta),
                0, 0, target_state_.v*std::cos(target_state_.theta),
                0, 0, 0;
        Eigen::Matrix<double, 3, 2> B;
        B<<std::cos(target_state_.theta), 0,
                std::sin(target_state_.theta), 0,
                0, 1;
        Eigen::Matrix<double, 3, 3> I = Eigen::Matrix<double, 3, 3>::Identity();
        A_ = (control_t * A + I);
        B_ = control_t * B;
        Q_ << 1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;
        R_ << 4, 0,
                    0, 4;
        // std::cout<<"******A******"<<std::endl;
        // std::cout<<A_<<std::endl;
        // std::cout<<"******B******"<<std::endl;
        // std::cout<<B_<<std::endl;
    }

    void setState(point cur_state, point target_state)
    {
        cur_state_ = cur_state;
        target_state_ = target_state;
    }

    void setControlRate(const double dt) {control_t = dt; }

    void setMaxIterCount(const int max_iter) { max_iter_count_ = max_iter; }

    void solve()
    {
        Eigen::Matrix<double, 3, 3> P;
        Eigen::Matrix<double, 3, 3> P_Next;
        P = Q_;
        int iter = 0;
        double eps = std::numeric_limits<double>::max();
        double eps_threshold = 1e-3;
        while (iter++ < max_iter_count_ && eps > eps_threshold)
        {
            P_Next = Q_ + A_.transpose() * P * A_ - A_.transpose() * P * B_* (R_ + B_.transpose()*P*B_).inverse()*B_.transpose()*P*A_;
            // eps = (P_Next - P).lpNorm<Eigen::Infinity>();
            eps = fabs((P_Next - P).maxCoeff());
            P = P_Next;
        }
        // std::cout<<"*****P****"<<std::endl;
        // std::cout<<P<<std::endl;
        printf("iter: %d, eps:%f.\r\n", iter, eps);

        K_ = -(R_ + B_.transpose()* P * B_).inverse()*B_.transpose()*P*A_;
        Eigen::Matrix<double, 3, 1> X;
        X << cur_state_.x - target_state_.x, cur_state_.y - target_state_.y, cur_state_.theta - target_state_.theta;
        u_ = K_ * X;
        ctrl.v = u_(0) + target_state_.v;
        ctrl.w = u_(1) + target_state_.w;
        // std::cout<<"*****X*****"<<std::endl;
        // std::cout<<X<<std::endl;
        // std::cout<<"*****K*****"<<std::endl;
        // std::cout<<K_<<std::endl;
        // std::cout<<"*****u*****"<<std::endl;
        // std::cout<<u_<<std::endl;
    }

    control getControlInput () const { return ctrl; }

private:
    Eigen::Matrix<double, 3, 3> A_;
    Eigen::Matrix<double, 3, 2> B_;
    Eigen::Matrix<double, 3, 3> Q_;
    Eigen::Matrix<double, 2, 2> R_;
    Eigen::Matrix<double, 2, 3> K_;
    Eigen::Matrix<double, 2, 1> u_;
    point cur_state_, target_state_;
    control ctrl;
    double control_t ;
    int max_iter_count_;
};

std::vector<point> getCirclePoint(point start_point, point end_point, double r)
{
    std::vector<point> res_points;
    for(double theta = start_point.theta; theta <= end_point.theta; theta = normalizeAngle(theta + M_PI/180))
    {
        point temp_point(start_point.x + r*std::cos(theta), start_point.y + r * std::sin(theta), normalizeAngle(theta + M_PI_2), 0, 0);
        if(res_points.empty() || temp_point != res_points.back()){
            res_points.emplace_back(temp_point);
            std::cout<<temp_point.x<<", "<<temp_point.y<<", "<<temp_point.theta<<std::endl;
        }
    }

    return res_points;
}

int getNearestId(const std::vector<point> &points, const point &point)
{
    int min_id = -1;
    double min_dis = std::numeric_limits<double>::max();
    int size = points.size();
    for(int i = 0; i < size; ++i){
        double dis = pow2(points[i].x - point.x) + pow2(points[i].y - point.y);
        if(dis < min_dis){
            min_dis = dis;
            min_id = i;
        }
    }
    return min_id;
}

std::vector<point> getLinePoints(point start_point, point end_point)
{
    std::vector<point> res_points;

    for(int i = 0; i < 100; ++i)
    {
        point temp_point(start_point.x + i * 0.01, 0, 0, 0, 0);
        res_points.emplace_back(temp_point);
        std::cout<<temp_point.x<<", "<<temp_point.y<<", "<<temp_point.theta<<std::endl;
    }

    return res_points;
}

int main()
{
    double dt = 0.05;
    point cur_state;
    std::shared_ptr<lqr> lqr_ptr = std::make_shared<lqr>();
    lqr_ptr->setControlRate(dt);
    lqr_ptr->setMaxIterCount(1000);

    std::vector<point> res_points = getCirclePoint(point(0, 10, -M_PI_2, 0, 0), point(1, 1, 0, 0, 0), 10);
    // std::vector<point> res_points = getLinePoints(point(0,0,0,0,0), point(0,0,0,0,0));
    int index = 0;
    while(index++ < 10 || true)
    {
        double dis = std::hypot(cur_state.x - res_points.back().x, cur_state.y - res_points.back().y);
        if(dis < 0.05){
            break;
        }
        int min_id = getNearestId(res_points, cur_state);
        min_id = min_id < res_points.size() - 1 ? min_id + 1 : min_id;
        point target_state = res_points[min_id];
        target_state.v = 0.5;
        target_state.w = normalizeAngle(target_state.theta - cur_state.theta);
        lqr_ptr->setState(cur_state, target_state);
        lqr_ptr->constructMatrix();
        lqr_ptr->solve();
        control ctrl = lqr_ptr->getControlInput();
        // cur_state.x = cur_state.x + ctrl.v/ctrl.w*(std::sin(cur_state.theta + ctrl.w * dt) - std::sin(cur_state.theta));
        // cur_state.y = cur_state.y + ctrl.v/ctrl.w*(std::cos(cur_state.theta + ctrl.w * dt) - std::cos(cur_state.theta));
        cur_state.x = cur_state.x + ctrl.v * dt * std::cos(cur_state.theta);
        cur_state.y = cur_state.y + ctrl.v * dt * std::sin(cur_state.theta);
        cur_state.theta = cur_state.theta + ctrl.w * dt;
        printf("index:%d, ref_x:%f, ref_y:%f, res_theta:%f, x:%f, y:%f, theta:%f, v:%f, w:%f.\r\n", 
            index, res_points[min_id].x, res_points[min_id].y, res_points[min_id].theta, cur_state.x, cur_state.y, cur_state.theta, ctrl.v, ctrl.w);
        if(ctrl.v < 1e-3 || ctrl.w < 1e-3)
        {
            break;
        }
        // std::cout<<"v:"<<ctrl.v<<", w:"<<ctrl.w<<std::endl;
    }
    
    
    return 0;
}
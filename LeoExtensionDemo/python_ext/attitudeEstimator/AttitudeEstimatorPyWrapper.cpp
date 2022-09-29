#include "AttitudeEstimatorPyWrapper.hpp"
#include "PyExt.hpp"
namespace gu
{
namespace attitudeEstimator
{
namespace py_ext
{

namespace py = boost::python;

WrappedMatrixReceiver::WrappedMatrixReceiver()
{
    add_listener(this);
    // listeners_list=py::list();
}
WrappedMatrixReceiver::~WrappedMatrixReceiver()
{
    for (auto it=pylistener_pointers.begin();it!=pylistener_pointers.end();++it)
    {
        py::decref(*it);
    }
}
void WrappedMatrixReceiver::add_python_listener(boost::python::object listener)
{
    PyObject* pt=listener.ptr();
    py::incref(pt);
    pylistener_pointers.push_back(pt);
}
void WrappedMatrixReceiver::handleUpdate(   const Matrix<double, 3, 7> &leds, 
                                            const Matrix<bool, 1, 7> &detected)
{
    gu::py_ext::Lock_PyGIL lock;
    try 
    {
        // std::cout<<"handling Update"<<std::endl;
        py::object py_leds = gu::py_ext::np_darray_from_eigen_matrix<3,7>(leds);
        py::object py_selectors =  gu::py_ext::np_barray_from_eigen_matrix<1,7>(detected);
        // std::cout<<"my list of listeners is "<<listeners_list<<std::endl;
        for (auto it=pylistener_pointers.begin();it!=pylistener_pointers.end();++it)
        {
            // std::cout<<"calling a callback function"<<std::endl;
            // py::object new_list=py::list();
            py::call<py::object>(*it, py_leds, py_selectors);
        }
        // std::cout<<"post calling that function"<<std::endl;
    }
    catch (const py::error_already_set &e){
        std::cout<< "Exception !"<<std::endl;
        try
        {
            PyErr_Print();
            py::object sys(py::handle<> (PyImport_ImportModule("sys")));
            py::object err = sys.attr("stderr");
            std::string err_text = py::extract<std::string>(err.attr("getvalue")());
            std::cerr<<err_text<<std::endl;
        } catch(...){
            std::cerr << "Failure to parse python error\n";
        }
    }
}

template<unsigned int N, unsigned int M>
Matrix<double,N,M> extract_dmatrix(py::object nparray)
{
    Matrix<double,N,M> mat;
    for (int i = 0;i<N;i++)
        for (int j=0;j<M;j++)
            mat(i,j)=py::extract<double>(nparray.slice(i,j));
    return mat;
}
template<unsigned int N, unsigned int M>
Matrix<bool,N,M> extract_bmatrix(py::object nparray)
{
    Matrix<bool,N,M> mat;
    for (int i = 0;i<N;i++)
        for (int j=0;j<M;j++)
            mat(i,j)=py::extract<char>(nparray.attr("__getitem__")(py::make_tuple(i,j)));
    return mat;
}
template<typename T, unsigned int N, unsigned int M>
Matrix<T,N,M> extract_matrix(py::object nparray)
{
    Matrix<T,N,M> mat;
    for (int i = 0;i<N;i++)
        for (int j=0;j<M;j++)
            mat(i,j)=py::extract<T>(nparray.attr("__getitem__")(py::make_tuple(i,j)))();
    return mat;
}

WrappedKalmanFilter::WrappedKalmanFilter(py::object camera_transform, 
    py::object constellation, py::object priorPointWeight, py::object priorMatrixWeight):
    gu::attitudeEstimator::HumePoseKalmanFilter(py::extract<double>(priorPointWeight)(), 
        py::extract<double>(priorMatrixWeight)(), Matrix3d::Identity(), Vector3d::Zero(), Vector3d::Zero() )
{
    Matrix<double,3,7> ledsInBody = extract_matrix<double,3,7>(constellation);
    Matrix3d cameraDCM=extract_matrix<double,3,3>(camera_transform);
    Vector3d ledCentroid = ledsInLEDFrame.rowwise().sum() * (1.0 / 7.0);
    setLEDPose(SpatialTransform<bodyFrame,ledFrame>());
    ledPose.setPosition(ledCentroid);
    ledsInLEDFrame=ledPose.inverse()*ledsInBody;
    cameraPose.setOrientation(cameraDCM);
    cameraPose.setPosition(Vector3d::Zero());
}
void WrappedKalmanFilter::measureLeds(py::object leds, py::object selector)
{
    Matrix<double,3,7> c_leds=extract_matrix<double,3,7>(leds);
    Matrix<bool,1,7> c_sel;
    for (int i=0;i<7;i++)
    {
        bool b = py::extract<bool>(selector[i])();
        c_sel(0,i)=b;
    }
    handleUpdate(c_leds,c_sel);
    // updateIMU(0.0,0.0,90.8,0.0,0.0,0.0);
}

py::object WrappedKalmanFilter::estimateTransform(void)
{
    py::object dcm = gu::py_ext::np_darray_from_eigen_matrix<3,3>(getTransform().getDCM());
    py::object vec =gu::py_ext::np_darray_from_eigen_matrix<3,1>(getTransform().getPosition());
    return py::make_tuple(dcm,vec);
}

void WrappedKalmanFilter::measureImu(py::object imu)
{
    updateIMU(
        py::extract<double>(imu[0])(),
        py::extract<double>(imu[1])(),
        py::extract<double>(imu[2])(),
        py::extract<double>(imu[3])(),
        py::extract<double>(imu[4])(),
        py::extract<double>(imu[5])());
}
py::object WrappedKalmanFilter::estimateLeds(void)
{
    SpatialTransform<worldFrame,bodyFrame> trans=getTransform();
    Matrix<double,3,7> leds=trans.getDCM()*ledsInLEDFrame;
    leds.colwise()+=trans.getPosition();
    py::object leds2 = gu::py_ext::np_darray_from_eigen_matrix<3,7>(leds);
    return leds2;
}

}
}
}
import vtk
import os

TEST_OUTPUT_FOLDER = os.path.join("tests", "test_output")

def test_file_writer_output():
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(os.path.join(TEST_OUTPUT_FOLDER, "test_file_writer.vtp"))
    reader.Update()

    pdata = reader.GetOutput()
    assert pdata.GetNumberOfCells() == 10
    assert pdata.GetNumberOfPoints() == 10
    assert pdata.GetPoint(5) == (5.0, 5.0, 0.0)
    assert (pdata.GetPointData().GetArray('Velocity').GetTuple(5)
            == (5.0, -5.0, 0.0))
    assert (pdata.GetPointData().GetArray('Pressure').GetValue(5)
            == 5.0)
    

def test_SPH_2D_output():
    for vtp_file in os.listdir(TEST_OUTPUT_FOLDER):
        file_path = os.path.join(TEST_OUTPUT_FOLDER, vtp_file)
        if os.path.isdir(file_path) or "test_file_writer" in vtp_file or not vtp_file.endswith("vtp"):
            continue
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(file_path)
        reader.Update()
        pdata = reader.GetOutput()
        assert pdata.GetNumberOfPoints() == 2182
        assert check_boundary(pdata.GetNumberOfPoints(), pdata.GetPointData().GetArray('Velocity'), pdata.GetPointData().GetArray('Boundary'))
        assert check_points(pdata.GetNumberOfPoints(), pdata, pdata.GetPointData().GetArray('Boundary'))
        assert check_velocity(pdata.GetNumberOfPoints(), pdata.GetPointData().GetArray('Velocity'))

def check_velocity(num_points, v_arr):
    for i in range(num_points):
        vel = v_arr.GetTuple(i)
        if abs(vel[0]) >= 20 or abs(vel[1]) >= 20:
            print("Velocity of fluid particles should be less than c0.")
            return False
    return True

def check_points(num_points, p_arr, b_arr):
    for i in range(num_points):
        pos = p_arr.GetPoint(i)
        if b_arr.GetValue(i) == 0:
            if pos[0] < 0 or pos[0] > 20 or pos[1] < 0 or pos[1] > 10:
                print("Fluid particles should be inside boundaries.")
                return False
    return True

def check_boundary(num_points, v_arr, b_arr):
    for i in range(num_points):
        if b_arr.GetValue(i) == 1:
            vel = v_arr.GetTuple(i)
            if abs(vel[0]) != 0 or abs(vel[1]) != 0:
                print("Velocity of boundary particles should be zero.")
                return False
    return True

clear all; clc;
nx = 8;
ny = 8;
nu = 4;
nd = 4;
tt = [1:1:nu]; td = [nu+1:1:nu+nd];
nlobj = nlmpc(nx,ny,'MV',tt,'MD',td);
%--------------------------------------------------------------------------
Ts = 0.05;
PredictionHorizon = 5;
nlobj.Ts = Ts;
nlobj.PredictionHorizon = PredictionHorizon;
nlobj.ControlHorizon = PredictionHorizon;
%---------------------------Plant Model-----------------------------------------------
nlobj.Model.StateFcn = @cds_model;
nlobj.Model.OutputFcn = @modelOutputFcn;
nlobj.Jacobian.OutputFcn = @JaccoModelOutputFcn;
%---------------------------Cost Fcn--------------------------------------
nlobj.Weights.OutputVariables = [10000 10000 10000 20000 1000 1000 1000 1];
nlobj.Weights.ManipulatedVariables = [1 1 1 0.1];
nlobj.Weights.ManipulatedVariablesRate = [0.1 0.1 0.1 0];
%---------------------------Constraints--------------------------------------

nlobj.Optimization.CustomIneqConFcn = @myIneqConFunction;
%--------------------------Optimization Specification------------------------------------------------
nlobj.Optimization.SolverOptions.Algorithm = 'sqp';
nlobj.Optimization.SolverOptions.MaxIterations = 400;
nlobj.Optimization.SolverOptions.StepTolerance = 1e-6;
nlobj.Optimization.SolverOptions.ConstraintTolerance = 1e-6;
nlobj.Optimization.SolverOptions.OptimalityTolerance = 1e-6;
%------------------------------------------------------------------
x0 = [0.45;0.7;0.5;0;0;0;0;0];
u0 = [0;0;0;0];
md0 = [0; 0; 0; 0];
validateFcns(nlobj,x0,u0,md0',[]);

Duration = 5;
Tsteps = Duration/Ts;
xHistory = x0';
Xref  = create_reference(Ts,Duration,PredictionHorizon);  
time = (0:Ts:(Duration+PredictionHorizon*Ts))';
ref.time  = time;
ref.signals.values = Xref;





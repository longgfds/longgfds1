using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using W_Entity;
using WZ_CommonLib;

namespace HMXS_MotionSys
{
    public class Algorithm
    {
        private static double M { get { return CommonMethod.flowPrmInfo.O_L_HorDistL_m; } }
        private static double N { get { return CommonMethod.flowPrmInfo.O_L_HorDistB_n; } }
        private static double DisRY { get { return CommonMethod.flowPrmInfo.R1_Y2N_Dist + CommonMethod.flowPrmInfo.Y2_Y1N_Dist + CommonMethod.flowPrmInfo.Y1_O_Dist; } }
        //侧板倾斜角度
        private static double UnderAngle { get { return CommonMethod.flowPrmInfo.UnderSideAngle; } }
        private static double UpAngle { get { return CommonMethod.flowPrmInfo.UpSideAngle; } }
        private static double OffUnderAngle { get { return CommonMethod.flowPrmInfo.OfflineSideAngle; } }

        //从R3旋转点到焦点距离
        private static double DisFocus { get { return CommonMethod.flowPrmInfo.R3_LaserV + CommonMethod.flowPrmInfo.LaserH; } }
        private static double ULDisFocus { get { return CommonMethod.flowPrmInfo.R3_LaserV + CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.ULZoomOff; } }
        private static double DRDisFocus { get { return CommonMethod.flowPrmInfo.R3_LaserV + CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.DRZoomOff; } }
        private static double OffULDisFocus { get { return CommonMethod.flowPrmInfo.R3_LaserV + CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.OfflineULZoomOff; } }
        private static double OffDRDisFocus { get { return CommonMethod.flowPrmInfo.R3_LaserV + CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.OfflineDRZoomOff; } }

        private static double unit = 0.5;
        public static int NumPoints { get { return Convert.ToInt32((360 + ExtraDegrees) / unit); } }
        private static int ExtraDegrees { get { return CommonMethod.flowPrmInfo.MarginAngle; } }

        public static double DownY1Off { get; set; }
        public static double UpY1Off { get; set; }

        private const short PT_SEGMENT_NORMAL = 0;
        private const short PT_SEGMENT_EVEN = 1;
        private const short PT_SEGMENT_STOP = 2;

        public static void CalculateDown(FourDimensionalCoordinate coordinate, ref WorkPieceInfo workP)
        {
            workP.SaveCentring = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            //其他保存赋初值
            workP.SaveDownStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SaveUpStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SavePartternPrFirstStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SavePartternPrSecondStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SavePartternPrThirdStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.OfflineSavePartternPrFirstStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.OfflineSavePartternPrSecondStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.OfflineSavePartternPrThirdStart = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SaveDownStart_B1 = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SaveUpStart_B1 = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SaveDownStart_B2 = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);
            workP.SaveUpStart_B2 = new FourDimensionalCoordinate(coordinate.R, coordinate.Y2, coordinate.Z, coordinate.Y1);


            FollowCal(ref workP);
            double xb = Math.Sqrt(Pow(N) + Pow(M));
            double a1 = RA(Math.Asin(N / xb));
            double a2 = RA(Math.Acos(N / xb));
            double by = DisRY + coordinate.Y1 + coordinate.Y2;
            double px = coordinate.R <= 90 ? -by * Math.Cos(AR(coordinate.R)) : by * Math.Cos(AR(180 - coordinate.R));
            double py = by * Math.Sin(AR(coordinate.R));

            //斜a  Rmb Rnc   ab夹角=90-asin(n/xb)
            double ang = 90 - a1;
            double c = Math.Abs(Math.Sqrt(Pow(xb) + Pow(by) - 2 * xb * by * Math.Cos(AR(ang))));
            //bc
            double angα = RA(Math.Acos((Pow(by) + Pow(c) - Pow(xb)) / (2 * by * c)));


            double fx, fy;
            if (coordinate.R + angα < 90)
            {
                fx = -c * Math.Cos(AR(coordinate.R + angα));
                fy = c * Math.Sin(AR(coordinate.R + angα));

            }
            else
            {
                fx = c * Math.Cos(AR(180 - coordinate.R - angα));
                fy = c * Math.Sin(AR(180 - coordinate.R - angα));
            }

            //执行半径 
            double clearR = workP.HubDiameter / 2.0 + workP.RingWidth / 2.0;
            double clearR_B1 = workP.HubDiameter / 2.0 + workP.RingWidth * 2 / 3.0 + CommonMethod.flowPrmInfo.Big_ClearOff1;
            double clearR_B2 = workP.HubDiameter / 2.0 + workP.RingWidth / 3.0 + CommonMethod.flowPrmInfo.Big_ClearOff2;

            //旋转偏移量  
            Tuple<double, double> downTuple;
            Tuple<double, double> downTuple_B1;
            Tuple<double, double> downTuple_B2;
            if (workP.ProductType == "线下侧板")
            {
                downTuple = DownCB_R3Rota_OffYZ(-OffUnderAngle);
                downTuple_B1 = DownCB_R3Rota_OffYZ(-CommonMethod.flowPrmInfo.Big_OfflineSideAngle1);
                downTuple_B2 = DownCB_R3Rota_OffYZ(-CommonMethod.flowPrmInfo.Big_OfflineSideAngle2);
            }
            else
            {
                downTuple = DownCB_R3Rota_OffYZ(-UnderAngle);
                downTuple_B1 = DownCB_R3Rota_OffYZ(-CommonMethod.flowPrmInfo.Big_UnderSideAngle1);
                downTuple_B2 = DownCB_R3Rota_OffYZ(-CommonMethod.flowPrmInfo.Big_UnderSideAngle2);
            }


            double pr = Math.Abs(Math.Sqrt(Pow(clearR + Math.Abs(downTuple.Item1) + N) + Pow(M)));
            double pr_B1 = Math.Abs(Math.Sqrt(Pow(clearR_B1 + Math.Abs(downTuple_B1.Item1) + N) + Pow(M)));
            double pr_B2 = Math.Abs(Math.Sqrt(Pow(clearR_B2 + Math.Abs(downTuple_B2.Item1) + N) + Pow(M)));

            //下侧板起始点极坐标Y
            double down_startxb = by + clearR + downTuple.Item1;
            double down_startpx = coordinate.R <= 90 ? -down_startxb * Math.Cos(AR(coordinate.R)) : down_startxb * Math.Cos(AR(180 - coordinate.R));
            double down_startpy = down_startxb * Math.Sin(AR(coordinate.R));
            double down_startz = workP.DownTTPM_CurZ + (CommonMethod.flowPrmInfo.LaserH - (workP.DownPM_Height - CommonMethod.flowPrmInfo.C_L_VerDist)) + downTuple.Item2;

            double down_ttxb = by - workP.HubDiameter / 2.0;
            double down_pmxb = by - workP.HubDiameter / 2.0 - workP.RingWidth / 2.0;

            double down_startxb_B1 = by + clearR_B1 + downTuple_B1.Item1;
            double down_startpx_B1 = coordinate.R <= 90 ? -down_startxb_B1 * Math.Cos(AR(coordinate.R)) : down_startxb_B1 * Math.Cos(AR(180 - coordinate.R));
            double down_startpy_B1 = down_startxb_B1 * Math.Sin(AR(coordinate.R));
            double down_startz_B1 = workP.DownTTPM_CurZ + (CommonMethod.flowPrmInfo.LaserH - (workP.DownPM_Height - CommonMethod.flowPrmInfo.C_L_VerDist)) + downTuple_B1.Item2;
            double down_startxb_B2 = by + clearR_B2 + downTuple_B2.Item1;
            double down_startpx_B2 = coordinate.R <= 90 ? -down_startxb_B2 * Math.Cos(AR(coordinate.R)) : down_startxb_B2 * Math.Cos(AR(180 - coordinate.R));
            double down_startpy_B2 = down_startxb_B2 * Math.Sin(AR(coordinate.R));
            double down_startz_B2 = workP.DownTTPM_CurZ + (CommonMethod.flowPrmInfo.LaserH - (workP.DownPM_Height - CommonMethod.flowPrmInfo.C_L_VerDist)) + downTuple_B2.Item2;

            workP.CenterCir.X = fx; workP.CenterCir.Y = fy;
            workP.DownCenter_PolarP.R = coordinate.R; workP.DownCenter_PolarP.Hy = by;
            workP.DownCenter_OrthP.X = px; workP.DownCenter_OrthP.Y = py;
            workP.DownStart_PolarP.R = coordinate.R; workP.DownStart_PolarP.Hy = down_startxb;
            workP.DownStart_OrthP.X = down_startpx; workP.DownStart_OrthP.Y = down_startpy;
            workP.DownStart_PolarP_B1.R = coordinate.R; workP.DownStart_PolarP_B1.Hy = down_startxb_B1;
            workP.DownStart_OrthP_B1.X = down_startpx_B1; workP.DownStart_OrthP_B1.Y = down_startpy_B1;
            workP.DownStart_PolarP_B2.R = coordinate.R; workP.DownStart_PolarP_B2.Hy = down_startxb_B2;
            workP.DownStart_OrthP_B2.X = down_startpx_B2; workP.DownStart_OrthP_B2.Y = down_startpy_B2;

            workP.DownTT_PolarP.R = coordinate.R; workP.DownTT_PolarP.Hy = down_ttxb;
            workP.DownPM_PolarP.R = coordinate.R; workP.DownPM_PolarP.Hy = down_pmxb;
            workP.DownY1Off = downTuple.Item1;
            workP.DownY1Off_B1 = downTuple_B1.Item1;
            workP.DownY1Off_B2 = downTuple_B2.Item1;
            workP.Start_Rmn_A = angα;

            workP.CBFr = clearR; workP.CBPr = pr;
            workP.CBFr_B1 = clearR_B1; workP.CBPr_B1 = pr_B1;
            workP.CBFr_B2 = clearR_B2; workP.CBPr_B2 = pr_B2;

            workP.DownWorkZ = down_startz;
            workP.DownWorkZ_B1 = down_startz_B1;
            workP.DownWorkZ_B2 = down_startz_B2;

            workP.SaveDownStart.Y1 = workP.SaveDownStart.Y1 + downTuple.Item1;
            workP.SaveDownStart.Y2 = workP.SaveDownStart.Y2 + clearR;
            workP.SaveDownStart_B1.Y1 = workP.SaveDownStart_B1.Y1 + downTuple_B1.Item1;
            workP.SaveDownStart_B1.Y2 = workP.SaveDownStart_B1.Y2 + clearR_B1;
            workP.SaveDownStart_B2.Y1 = workP.SaveDownStart_B2.Y1 + downTuple_B2.Item1;
            workP.SaveDownStart_B2.Y2 = workP.SaveDownStart_B2.Y2 + clearR_B2;

            workP.DownLightWidth =
                workP.DownLightWidth_B1 =
                workP.DownLightWidth_B2 = (workP.RingWidth + CommonMethod.flowPrmInfo.CBLightOff) > CommonMethod.hZCardInfo.LightWidth ? CommonMethod.hZCardInfo.LightWidth : (workP.RingWidth + CommonMethod.flowPrmInfo.CBLightOff);

            workP.HWKLightWidth = (workP.TreadWidth + CommonMethod.flowPrmInfo.HWKLightOff) > CommonMethod.hZCardInfo.LightWidth ? CommonMethod.hZCardInfo.LightWidth : (workP.TreadWidth + CommonMethod.flowPrmInfo.HWKLightOff);
        }


        public static void CalculateUp(ref WorkPieceInfo workP)
        {
            Tuple<double, double> upTuple = UpCB_R3Rota_OffYZ(UpAngle);
            Tuple<double, double> upTuple_B1 = UpCB_R3Rota_OffYZ(CommonMethod.flowPrmInfo.Big_UpSideAngle1);
            Tuple<double, double> upTuple_B2 = UpCB_R3Rota_OffYZ(CommonMethod.flowPrmInfo.Big_UpSideAngle2);

            //上侧板起始点极坐标 直角坐标
            double up_centerxb = workP.DownCenter_PolarP.Hy - N * 2;
            double up_startxb = up_centerxb + workP.CBFr + upTuple.Item1;
            double up_startpx = workP.DownCenter_PolarP.R <= 90 ? -up_startxb * Math.Cos(AR(workP.DownCenter_PolarP.R)) : up_startxb * Math.Cos(AR(180 - workP.DownCenter_PolarP.R));
            double up_startpy = up_startxb * Math.Sin(AR(workP.DownCenter_PolarP.R));
            double up_startz = workP.UpTTPM_CurZ + ((workP.UpTT_Height + workP.UpPM_Height) / 2.0 - CommonMethod.flowPrmInfo.C_L_VerDist - CommonMethod.flowPrmInfo.LaserH) + upTuple.Item2;

            //上侧凸台坐标及侧板平面坐标(激光头垂直时的坐标  最终要用到上侧板的焦距和花纹块的中心计算)
            double up_ttxb = up_centerxb + workP.HubDiameter / 2.0;
            double up_pmxb = up_centerxb + workP.HubDiameter / 2.0 + workP.RingWidth / 2.0;

            double up_startxb_B1 = up_centerxb + workP.CBFr_B1 + upTuple_B1.Item1;
            double up_startpx_B1 = workP.DownCenter_PolarP.R <= 90 ? -up_startxb_B1 * Math.Cos(AR(workP.DownCenter_PolarP.R)) : up_startxb_B1 * Math.Cos(AR(180 - workP.DownCenter_PolarP.R));
            double up_startpy_B1 = up_startxb_B1 * Math.Sin(AR(workP.DownCenter_PolarP.R));
            double up_startz_B1 = workP.UpTTPM_CurZ - ((workP.UpTT_Height + workP.UpPM_Height) / 2.0 - CommonMethod.flowPrmInfo.C_L_VerDist - CommonMethod.flowPrmInfo.LaserH) + upTuple_B1.Item2;
            double up_startxb_B2 = up_centerxb + workP.CBFr_B2 + upTuple_B2.Item1;
            double up_startpx_B2 = workP.DownCenter_PolarP.R <= 90 ? -up_startxb_B2 * Math.Cos(AR(workP.DownCenter_PolarP.R)) : up_startxb_B2 * Math.Cos(AR(180 - workP.DownCenter_PolarP.R));
            double up_startpy_B2 = up_startxb_B2 * Math.Sin(AR(workP.DownCenter_PolarP.R));
            double up_startz_B2 = workP.UpTTPM_CurZ - ((workP.UpTT_Height + workP.UpPM_Height) / 2.0 - CommonMethod.flowPrmInfo.C_L_VerDist - CommonMethod.flowPrmInfo.LaserH) + upTuple_B2.Item2;

            workP.UpCenter_PolarP.R = workP.DownCenter_PolarP.R; workP.UpCenter_PolarP.Hy = up_centerxb;
            workP.UpStart_PolarP.R = workP.DownCenter_PolarP.R; workP.UpStart_PolarP.Hy = up_startxb;
            workP.UpStart_OrthP.X = up_startpx; workP.UpStart_OrthP.Y = up_startpy;
            workP.UpStart_PolarP_B1.R = workP.DownCenter_PolarP.R; workP.UpStart_PolarP_B1.Hy = up_startxb_B1;
            workP.UpStart_OrthP_B1.X = up_startpx_B1; workP.UpStart_OrthP_B1.Y = up_startpy_B1;
            workP.UpStart_PolarP_B2.R = workP.DownCenter_PolarP.R; workP.UpStart_PolarP_B2.Hy = up_startxb_B2;
            workP.UpStart_OrthP_B2.X = up_startpx_B2; workP.UpStart_OrthP_B2.Y = up_startpy_B2;

            workP.UpTT_PolarP.R = workP.DownCenter_PolarP.R; workP.UpTT_PolarP.Hy = up_ttxb;
            workP.UpPM_PolarP.R = workP.DownCenter_PolarP.R; workP.UpPM_PolarP.Hy = up_pmxb;
            workP.UpY1Off = upTuple.Item1;
            workP.UpY1Off_B1 = upTuple_B1.Item1;
            workP.UpY1Off_B2 = upTuple_B2.Item1;
            workP.UpWorkZ = up_startz;
            workP.UpWorkZ_B1 = up_startz_B1;
            workP.UpWorkZ_B2 = up_startz_B2;

            workP.SaveUpStart.Y1 = workP.SaveUpStart.Y1 + upTuple.Item1;
            workP.SaveUpStart.Y2 = workP.SaveUpStart.Y2 + workP.CBFr - 2 * N;
            workP.SaveUpStart_B1.Y1 = workP.SaveUpStart_B1.Y1 + upTuple_B1.Item1;
            workP.SaveUpStart_B1.Y2 = workP.SaveUpStart_B1.Y2 + workP.CBFr_B1 - 2 * N;
            workP.SaveUpStart_B2.Y1 = workP.SaveUpStart_B2.Y1 + upTuple_B2.Item1;
            workP.SaveUpStart_B2.Y2 = workP.SaveUpStart_B2.Y2 + workP.CBFr_B2 - 2 * N;



        }



        public static void CalculateTrajectory_DownCB(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownStart_OrthP.X, workP.DownStart_OrthP.Y, workP.SaveDownStart.Y2, workP.DownStart_PolarP.Hy);
            //Tuple<OrthogonalCoor, PolarCoor, double> coordinate = Tuple.Create<OrthogonalCoor, PolarCoor, double>(new OrthogonalCoor(workP.CenterCir.X, workP.CenterCir.Y),new PolarCoor(workP.DownCenter_PolarP.R, workP.DownCenter_PolarP.Hy),0);
            //workP.DownSideCoordinates.Add(coordinate);

            workP.DownSideCoordinates = workP.DownSideCoordinates.GetRange(workP.DownSideCoordinates.Count / 2, workP.DownSideCoordinates.Count - workP.DownSideCoordinates.Count / 2 - 1)
                                        .Concat(workP.DownSideCoordinates.GetRange(0, workP.DownSideCoordinates.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                        .ToList();
        }

        public static void CalculateTrajectory_DownCBL(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinatesL = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownLStart_OrthP.X, workP.DownLStart_OrthP.Y, workP.SaveDownLStart.Y2, workP.DownLStart_PolarP.Hy);
            //workP.DownSideCoordinatesL = workP.DownSideCoordinatesL.GetRange(workP.DownSideCoordinatesL.Count / 2, workP.DownSideCoordinatesL.Count - workP.DownSideCoordinatesL.Count / 2 - 1)
            //                            .Concat(workP.DownSideCoordinatesL.GetRange(0, workP.DownSideCoordinatesL.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
            //                            .ToList();
        }
        public static void CalculateTrajectory_DownCBR(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinatesR = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownRStart_OrthP.X, workP.DownRStart_OrthP.Y, workP.SaveDownRStart.Y2, workP.DownRStart_PolarP.Hy);
            //workP.DownSideCoordinatesR = workP.DownSideCoordinatesR.GetRange(workP.DownSideCoordinatesR.Count / 2, workP.DownSideCoordinatesR.Count - workP.DownSideCoordinatesR.Count / 2 - 1)
            //                            .Concat(workP.DownSideCoordinatesR.GetRange(0, workP.DownSideCoordinatesR.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
            //                            .ToList();
        }

        public static void CalculateTrajectory_DownCB_B1(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates_B1 = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownStart_OrthP_B1.X, workP.DownStart_OrthP_B1.Y, workP.SaveDownStart_B1.Y2, workP.DownStart_PolarP_B1.Hy);
            workP.DownSideCoordinates_B1 = workP.DownSideCoordinates_B1.GetRange(workP.DownSideCoordinates_B1.Count / 2, workP.DownSideCoordinates_B1.Count - workP.DownSideCoordinates_B1.Count / 2 - 1)
                                         .Concat(workP.DownSideCoordinates_B1.GetRange(0, workP.DownSideCoordinates_B1.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                         .ToList();
        }
        public static void CalculateTrajectory_DownCB_B2(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates_B2 = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownStart_OrthP_B2.X, workP.DownStart_OrthP_B2.Y, workP.SaveDownStart_B2.Y2, workP.DownStart_PolarP_B2.Hy);
            workP.DownSideCoordinates_B2 = workP.DownSideCoordinates_B2.GetRange(workP.DownSideCoordinates_B2.Count / 2, workP.DownSideCoordinates_B2.Count - workP.DownSideCoordinates_B2.Count / 2 - 1)
                                        .Concat(workP.DownSideCoordinates_B2.GetRange(0, workP.DownSideCoordinates_B2.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                        .ToList();
        }
        public static void CalculateTrajectory_DownCB_LB1(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates_LB1 = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownLStart_OrthP_B1.X, workP.DownLStart_OrthP_B1.Y, workP.SaveDownLStart_B1.Y2, workP.DownLStart_PolarP_B1.Hy);
            workP.DownSideCoordinates_LB1 = workP.DownSideCoordinates_B1.GetRange(workP.DownSideCoordinates_LB1.Count / 2, workP.DownSideCoordinates_LB1.Count - workP.DownSideCoordinates_LB1.Count / 2 - 1)
                                         .Concat(workP.DownSideCoordinates_LB1.GetRange(0, workP.DownSideCoordinates_LB1.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                         .ToList();
        }
        public static void CalculateTrajectory_DownCB_LB2(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates_LB2 = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownLStart_OrthP_B2.X, workP.DownLStart_OrthP_B2.Y, workP.SaveDownLStart_B2.Y2, workP.DownLStart_PolarP_B2.Hy);
            workP.DownSideCoordinates_LB2 = workP.DownSideCoordinates_LB2.GetRange(workP.DownSideCoordinates_LB2.Count / 2, workP.DownSideCoordinates_LB2.Count - workP.DownSideCoordinates_LB2.Count / 2 - 1)
                                         .Concat(workP.DownSideCoordinates_LB2.GetRange(0, workP.DownSideCoordinates_LB2.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                         .ToList();
        }
        public static void CalculateTrajectory_DownCB_RB1(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates_RB1 = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownRStart_OrthP_B1.X, workP.DownRStart_OrthP_B1.Y, workP.SaveDownRStart_B1.Y2, workP.DownRStart_PolarP_B1.Hy);
            workP.DownSideCoordinates_RB1 = workP.DownSideCoordinates_RB1.GetRange(workP.DownSideCoordinates_RB1.Count / 2, workP.DownSideCoordinates_RB1.Count - workP.DownSideCoordinates_RB1.Count / 2 - 1)
                                         .Concat(workP.DownSideCoordinates_RB1.GetRange(0, workP.DownSideCoordinates_RB1.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                         .ToList();
        }
        public static void CalculateTrajectory_DownCB_RB2(ref WorkPieceInfo workP)
        {
            workP.DownSideCoordinates_RB2 = GetSideTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.DownRStart_OrthP_B2.X, workP.DownRStart_OrthP_B2.Y, workP.SaveDownRStart_B2.Y2, workP.DownRStart_PolarP_B2.Hy);
            workP.DownSideCoordinates_RB2 = workP.DownSideCoordinates_RB2.GetRange(workP.DownSideCoordinates_RB2.Count / 2, workP.DownSideCoordinates_RB2.Count - workP.DownSideCoordinates_RB2.Count / 2 - 1)
                                         .Concat(workP.DownSideCoordinates_RB2.GetRange(0, workP.DownSideCoordinates_RB2.Count / 2 + 1 + CommonMethod.flowPrmInfo.MarginAngle * 2))
                                         .ToList();
        }



        public static void CalculateTrajectory_UpCB(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpStart_OrthP.X, workP.UpStart_OrthP.Y, workP.SaveUpStart.Y2, workP.UpStart_PolarP.Hy);

        }

        public static void CalculateTrajectory_UpCB_B1(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates_B1 = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpStart_OrthP_B1.X, workP.UpStart_OrthP_B1.Y, workP.SaveUpStart_B1.Y2, workP.UpStart_PolarP_B1.Hy);

        }
        public static void CalculateTrajectory_UpCB_B2(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates_B2 = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpStart_OrthP_B2.X, workP.UpStart_OrthP_B2.Y, workP.SaveUpStart_B2.Y2, workP.UpStart_PolarP_B2.Hy);

        }

        public static void CalculateTrajectory_UpCBL(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinatesL = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpLStart_OrthP.X, workP.UpLStart_OrthP.Y, workP.SaveUpLStart.Y2, workP.UpLStart_PolarP.Hy);

        }
        public static void CalculateTrajectory_UpCBR(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinatesR = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpRStart_OrthP.X, workP.UpRStart_OrthP.Y, workP.SaveUpRStart.Y2, workP.UpRStart_PolarP.Hy);

        }
        public static void CalculateTrajectory_UpCB_LB1(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates_LB1 = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpLStart_OrthP_B1.X, workP.UpLStart_OrthP_B1.Y, workP.SaveUpLStart_B1.Y2, workP.UpLStart_PolarP_B1.Hy);

        }
        public static void CalculateTrajectory_UpCB_RB1(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates_RB1 = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpRStart_OrthP_B1.X, workP.UpRStart_OrthP_B1.Y, workP.SaveUpRStart_B1.Y2, workP.UpRStart_PolarP_B1.Hy);

        }
        public static void CalculateTrajectory_UpCB_LB2(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates_LB2 = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpLStart_OrthP_B2.X, workP.UpLStart_OrthP_B2.Y, workP.SaveUpLStart_B2.Y2, workP.UpLStart_PolarP_B2.Hy);

        }
        public static void CalculateTrajectory_UpCB_RB2(ref WorkPieceInfo workP)
        {
            workP.UpSideCoordinates_RB2 = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpRStart_OrthP_B2.X, workP.UpRStart_OrthP_B2.Y, workP.SaveUpRStart_B2.Y2, workP.UpRStart_PolarP_B2.Hy);

        }


        public static void CalculatePattern(ref WorkPieceInfo workP)
        {           
            workP.PartternFr = workP.TreadDiameter / 2.0;
            double pattern_centerxb = workP.UpCenter_PolarP.Hy + CommonMethod.flowPrmInfo.R3_LaserV + N;

            double LHeighe = (workP.UpTT_Height + workP.UpPM_Height) / 2.0 - CommonMethod.flowPrmInfo.C_L_VerDist;
            double LHorH = LHeighe + CommonMethod.flowPrmInfo.R3_LaserV - N;
            double LHor_CH = LHorH - workP.TreadWidth / 2.0;

            double pattern_centerZ = workP.UpTTPM_CurZ + LHor_CH;//向上为正 

            double LHeighe2 = workP.UpTT_Height - CommonMethod.flowPrmInfo.C_L_VerDist;
            double LHorH2 = LHeighe2 + CommonMethod.flowPrmInfo.R3_LaserV - N;
            double pattern_ttZ = workP.UpTTPM_CurZ + LHorH2 - CommonMethod.flowPrmInfo.NCL_Dis;

            //水平焦距时hy
            double pattern_zoomxb = pattern_centerxb + (CommonMethod.flowPrmInfo.LaserH - workP.TreadDiameter / 2.0);//小远离  大靠近
            double upattern_zoomxb = pattern_centerxb + (CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.ULZoomOff - workP.TreadDiameter / 2.0);
            double dpattern_zoomxb = pattern_centerxb + (CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.DRZoomOff - workP.TreadDiameter / 2.0);
            //花纹块左上坐标计算
            Tuple<double, double> upPatternTuple = Pattern_R3Rota_OffYZ(-workP.UA, 1);

            double upAngle_centerZ = 0;

            if (workP.TreadWidth > CommonMethod.flowPrmInfo.FourClearWidth)
                upAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.PULZOff2 : pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.PULZOff_N2;
            else
                upAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.PULZOff : pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.PULZOff_N;

            double upAngle_zoomxb = upattern_zoomxb + upPatternTuple.Item2;

            double ra = workP.DownCenter_PolarP.R;
            workP.UpHorCenter_PolarP.R = ra; workP.UpHorCenter_PolarP.Hy = pattern_centerxb;
            workP.UpHorZoom_PolarP.R = ra; workP.UpHorZoom_PolarP.Hy = pattern_zoomxb;
            workP.Pattern_CenterZ = pattern_centerZ;
            workP.UpTT_Z = pattern_ttZ;

            workP.UpPattern_PolarP.R = ra; workP.UpPattern_PolarP.Hy = upAngle_zoomxb;
            workP.UpPattern_OrthogonalP.X = Polar_Orthogonal(workP.UpPattern_PolarP).Item1;
            workP.UpPattern_OrthogonalP.Y = Polar_Orthogonal(workP.UpPattern_PolarP).Item2;
            workP.UpPattern_Z = upAngle_centerZ;

            //左夹角实现计算
            LAngleRealize(ref workP);


            //花纹块右下坐标计算
            Tuple<double, double> downPatternTuple = Pattern_R3Rota_OffYZ(workP.DA, 2);

            double downAngle_centerZ = 0;
            if (workP.TreadWidth > CommonMethod.flowPrmInfo.FourClearWidth)
                downAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.PDRZOff2 : pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.PDRZOff_N2;
            else
                downAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.PDRZOff : pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.PDRZOff_N;
            double downAngle_zoomxb = dpattern_zoomxb + downPatternTuple.Item2;

            workP.DownPattern_PolarP.R = ra; workP.DownPattern_PolarP.Hy = downAngle_zoomxb;
            workP.DownPattern_OrthogonalP.X = Polar_Orthogonal(workP.DownPattern_PolarP).Item1;
            workP.DownPattern_OrthogonalP.Y = Polar_Orthogonal(workP.DownPattern_PolarP).Item2;
            workP.DownPattern_Z = downAngle_centerZ;

            //右夹角实现计算
            RAngleRealize(ref workP);


            //特殊洗坐标计算  
            if (workP.HubDiameter <= CommonMethod.flowPrmInfo.HubLimit * 25.4)
            {
                workP.IsSmall = true;
                workP.SPSmall_Pattern_PolarP.R = ra; workP.SPSmall_Pattern_PolarP.Hy = pattern_zoomxb;
                workP.SPSmallPattern_OrthogonalP.X = Polar_Orthogonal(workP.SPSmall_Pattern_PolarP).Item1;
                workP.SPSmallPattern_OrthogonalP.Y = Polar_Orthogonal(workP.SPSmall_Pattern_PolarP).Item2;
                workP.SmallPattern_Z = pattern_ttZ + CommonMethod.flowPrmInfo.PSPZOff;
            }
            else
            {
                workP.SPBig_Pattern_PolarP.R = ra; workP.SPBig_Pattern_PolarP.Hy = pattern_centerxb;
                workP.SPBigPattern_OrthogonalP.X = Polar_Orthogonal(workP.SPBig_Pattern_PolarP).Item1;
                workP.SPBigPattern_OrthogonalP.Y = Polar_Orthogonal(workP.SPBig_Pattern_PolarP).Item2;
                workP.BigPattern_Z = pattern_ttZ + (workP.UpPM_Height - workP.UpTT_Height) + CommonMethod.flowPrmInfo.PSPZOff;
                //中心时N外端点就是往里走一个激光头外焦距
                workP.SPN_PolarP.R = ra; workP.SPN_PolarP.Hy = pattern_centerxb - CommonMethod.flowPrmInfo.LaserH;
                workP.SPN_OrthogonalP.X = Polar_Orthogonal(workP.SPN_PolarP).Item1;
                workP.SPN_OrthogonalP.Y = Polar_Orthogonal(workP.SPN_PolarP).Item2;

                if (!OAngleCalculate(ref workP).IsSuccess)//无解则使用小的
                {
                    workP.IsSmall = true;
                    workP.SPSmall_Pattern_PolarP.R = ra; workP.SPSmall_Pattern_PolarP.Hy = pattern_zoomxb;
                    workP.SPSmallPattern_OrthogonalP.X = Polar_Orthogonal(workP.SPSmall_Pattern_PolarP).Item1;
                    workP.SPSmallPattern_OrthogonalP.Y = Polar_Orthogonal(workP.SPSmall_Pattern_PolarP).Item2;
                    workP.SmallPattern_Z = pattern_ttZ;
                }
                else
                    workP.IsSmall = false;
            }


            //先到上中心Y2 再到水平中心Y2 再到水平焦距Y2 再到左上Y2
            workP.SavePartternPrFirstStart.R = workP.UpLeftPattern_PolarP.R;
            workP.SavePartternPrFirstStart.Y2 = workP.SavePartternPrFirstStart.Y2 + (workP.UpLeftPattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);

            workP.SavePartternPrSecondStart.R = workP.DownRightPattern_PolarP.R;
            workP.SavePartternPrSecondStart.Y2 = workP.SavePartternPrSecondStart.Y2 + (workP.DownRightPattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);

            if (workP.IsSmall)
            {
                workP.SavePartternPrThirdStart.Y2 = workP.SavePartternPrThirdStart.Y2 + (workP.SPSmall_Pattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);
            }
            else
            {
                workP.SavePartternPrThirdStart.Y2 = workP.SavePartternPrThirdStart.Y2 + (workP.SPBig_Pattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);
            }


        }



        public static void CalculateTrajectory_PatternUL(ref WorkPieceInfo workP)
        {
            double radius = Math.Sqrt(Pow(workP.CenterCir.X - workP.UpLeftPattern_OrthogonalP.X) + Pow(workP.CenterCir.Y - workP.UpLeftPattern_OrthogonalP.Y));
            workP.PartternPrFirst = radius;
            workP.ULPartternCoordinates = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.UpLeftPattern_OrthogonalP.X, workP.UpLeftPattern_OrthogonalP.Y, workP.SavePartternPrFirstStart.Y2, workP.UpLeftPattern_PolarP.Hy);

        }


        public static void CalculateTrajectory_PatternDR(ref WorkPieceInfo workP)
        {
            double radius = Math.Sqrt(Pow(workP.CenterCir.X - workP.DownRightPattern_OrthogonalP.X) + Pow(workP.CenterCir.Y - workP.DownRightPattern_OrthogonalP.Y));
            workP.PartternPrSecond = radius;
            workP.DRPartternCoordinates = GetTrajectoryF(workP.CenterCir.X, workP.CenterCir.Y, workP.DownRightPattern_OrthogonalP.X, workP.DownRightPattern_OrthogonalP.Y, workP.SavePartternPrSecondStart.Y2, workP.DownRightPattern_PolarP.Hy);

        }


        public static void CalculateTrajectory_PatternSp(ref WorkPieceInfo workP)
        {

            List<Tuple<OrthogonalCoor, PolarCoor, double>> coordinates = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();

            double x1, y1, xb;

            if (workP.IsSmall)
            {
                x1 = workP.SPSmallPattern_OrthogonalP.X;
                y1 = workP.SPSmallPattern_OrthogonalP.Y;
                xb = workP.SPSmall_Pattern_PolarP.Hy;
            }
            else
            {
                x1 = workP.SPBigPattern_OrthogonalP.X;
                y1 = workP.SPBigPattern_OrthogonalP.Y;
                xb = workP.SPBig_Pattern_PolarP.Hy;
            }
            double radius = Math.Sqrt(Pow(workP.CenterCir.X - x1) + Pow(workP.CenterCir.Y - y1));
            workP.PartternPrThird = radius;
            workP.SPPartternCoordinates = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, x1, y1, workP.SavePartternPrThirdStart.Y2, xb);


        }



        //OfflineParttern
        #region OfflineParttern
        public static void OfflineCalculatePattern(ref WorkPieceInfo workP)
        {
            workP.OfflinePartternFr = workP.TreadDiameter / 2.0;
            double pattern_centerxb = workP.DownCenter_PolarP.Hy + CommonMethod.flowPrmInfo.R3_LaserV - N;

            //计算花纹块中心Z
            double LHeighe = workP.DownPM_Height - CommonMethod.flowPrmInfo.C_L_VerDist;
            double LHorH = LHeighe + CommonMethod.flowPrmInfo.R3_LaserV + N;
            double LHor_CH = LHorH - workP.PatternH / 2.0;

            double pattern_centerZ = workP.DownTTPM_CurZ - LHor_CH;

            double LHeighe2 = workP.DownPM_Height - CommonMethod.flowPrmInfo.C_L_VerDist;
            double LHorH2 = LHeighe2 + CommonMethod.flowPrmInfo.R3_LaserV + N;
            double pattern_ttZ = workP.DownTTPM_CurZ - LHorH2 + CommonMethod.flowPrmInfo.CL_Dis;

            //水平焦距时hy
            double pattern_zoomxb = pattern_centerxb + (CommonMethod.flowPrmInfo.LaserH - workP.TreadDiameter / 2.0);
            double upattern_zoomxb = pattern_centerxb + (CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.OfflineULZoomOff - workP.TreadDiameter / 2.0);
            double dpattern_zoomxb = pattern_centerxb + (CommonMethod.flowPrmInfo.LaserH + CommonMethod.flowPrmInfo.OfflineDRZoomOff - workP.TreadDiameter / 2.0);//小远离  大靠近

            //花纹块左上坐标计算
            Tuple<double, double> upPatternTuple = Pattern_R3Rota_OffYZ(-CommonMethod.workPieceInfo.UA, 1);

            double upAngle_centerZ = 0;
            if (workP.TreadWidth > CommonMethod.flowPrmInfo.FourClearWidth)
                upAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineULZOff2 : pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineULZOff_N2;
            else
                upAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineULZOff : pattern_centerZ + upPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineULZOff_N;
            double upAngle_zoomxb = upattern_zoomxb + upPatternTuple.Item2;

            double ra = workP.DownCenter_PolarP.R;
            workP.OfflineHorCenter_PolarP.R = ra; workP.OfflineHorCenter_PolarP.Hy = pattern_centerxb;
            workP.OfflineHorZoom_PolarP.R = ra; workP.OfflineHorZoom_PolarP.Hy = pattern_zoomxb;
            workP.OfflinePattern_CenterZ = pattern_centerZ;
            workP.OfflineTT_Z = pattern_ttZ;

            workP.OfflineUpPattern_PolarP.R = ra; workP.OfflineUpPattern_PolarP.Hy = upAngle_zoomxb;
            workP.OfflineUpPattern_OrthogonalP.X = Polar_Orthogonal(workP.OfflineUpPattern_PolarP).Item1;
            workP.OfflineUpPattern_OrthogonalP.Y = Polar_Orthogonal(workP.OfflineUpPattern_PolarP).Item2;
            workP.OfflineUpPattern_Z = upAngle_centerZ;

            //左夹角实现计算
            OfflineLAngleRealize(ref workP);


            //花纹块右下坐标计算
            Tuple<double, double> downPatternTuple = Pattern_R3Rota_OffYZ(CommonMethod.workPieceInfo.DA, 2);

            double downAngle_centerZ = 0;
            if (workP.TreadWidth > CommonMethod.flowPrmInfo.FourClearWidth)
                downAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineDRZOff2 : pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineDRZOff_N2;
            else
                downAngle_centerZ = (workP.UA >= 0 && workP.DA >= 0) ? pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineDRZOff : pattern_centerZ + downPatternTuple.Item1 + CommonMethod.flowPrmInfo.OfflineDRZOff_N;
            double downAngle_zoomxb = dpattern_zoomxb + downPatternTuple.Item2;

            workP.OfflineDownPattern_PolarP.R = ra; workP.OfflineDownPattern_PolarP.Hy = downAngle_zoomxb;
            workP.OfflineDownPattern_OrthogonalP.X = Polar_Orthogonal(workP.OfflineDownPattern_PolarP).Item1;
            workP.OfflineDownPattern_OrthogonalP.Y = Polar_Orthogonal(workP.OfflineDownPattern_PolarP).Item2;
            workP.OfflineDownPattern_Z = downAngle_centerZ;

            //右夹角实现计算
            OfflineRAngleRealize(ref workP);


            //特殊洗坐标计算  
            if (true)//workP.HubDiameter <= CommonMethod.flowPrmInfo.HubLimit * 25.4
            {
                workP.OfflineIsSmall = true;
                workP.OfflineSPSmall_Pattern_PolarP.R = ra; workP.OfflineSPSmall_Pattern_PolarP.Hy = pattern_zoomxb;
                workP.OfflineSPSmallPattern_OrthogonalP.X = Polar_Orthogonal(workP.OfflineSPSmall_Pattern_PolarP).Item1;
                workP.OfflineSPSmallPattern_OrthogonalP.Y = Polar_Orthogonal(workP.OfflineSPSmall_Pattern_PolarP).Item2;
                workP.OfflineSmallPattern_Z = pattern_ttZ + CommonMethod.flowPrmInfo.OfflineSPZOff;
            }
            else
            {
                workP.OfflineSPBig_Pattern_PolarP.R = ra; workP.OfflineSPBig_Pattern_PolarP.Hy = pattern_centerxb;
                workP.OfflineSPBigPattern_OrthogonalP.X = Polar_Orthogonal(workP.OfflineSPBig_Pattern_PolarP).Item1;
                workP.OfflineSPBigPattern_OrthogonalP.Y = Polar_Orthogonal(workP.OfflineSPBig_Pattern_PolarP).Item2;
                workP.OfflineBigPattern_Z = pattern_ttZ - (workP.DownPM_Height - workP.DownTT_Height);
                //中心时N外端点就是往里走一个激光头外焦距
                workP.OfflineSPN_PolarP.R = ra; workP.OfflineSPN_PolarP.Hy = pattern_centerxb - CommonMethod.flowPrmInfo.LaserH;
                workP.OfflineSPN_OrthogonalP.X = Polar_Orthogonal(workP.OfflineSPN_PolarP).Item1;
                workP.OfflineSPN_OrthogonalP.Y = Polar_Orthogonal(workP.OfflineSPN_PolarP).Item2;

                if (!OAngleCalculate(ref workP).IsSuccess)//无解则使用小的
                {
                    workP.OfflineIsSmall = true;
                    workP.OfflineSPSmall_Pattern_PolarP.R = ra; workP.OfflineSPSmall_Pattern_PolarP.Hy = pattern_zoomxb;
                    workP.OfflineSPSmallPattern_OrthogonalP.X = Polar_Orthogonal(workP.OfflineSPSmall_Pattern_PolarP).Item1;
                    workP.OfflineSPSmallPattern_OrthogonalP.Y = Polar_Orthogonal(workP.OfflineSPSmall_Pattern_PolarP).Item2;
                    workP.OfflineSmallPattern_Z = pattern_ttZ;
                }
                else
                    workP.OfflineIsSmall = false;
            }


            //先到上中心Y2 再到水平中心Y2 再到水平焦距Y2 再到左上Y2
            workP.OfflineSavePartternPrFirstStart.R = workP.OfflineUpLeftPattern_PolarP.R;
            workP.OfflineSavePartternPrFirstStart.Y2 = workP.OfflineSavePartternPrFirstStart.Y2 + (workP.OfflineUpLeftPattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);

            workP.OfflineSavePartternPrSecondStart.R = workP.OfflineDownRightPattern_PolarP.R;
            workP.OfflineSavePartternPrSecondStart.Y2 = workP.OfflineSavePartternPrSecondStart.Y2 + (workP.OfflineDownRightPattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);

            if (workP.OfflineIsSmall)
            {
                workP.OfflineSavePartternPrThirdStart.Y2 = workP.OfflineSavePartternPrThirdStart.Y2 + (workP.OfflineSPSmall_Pattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);
            }
            else
            {
                workP.OfflineSavePartternPrThirdStart.Y2 = workP.OfflineSavePartternPrThirdStart.Y2 + (workP.OfflineSPBig_Pattern_PolarP.Hy - workP.DownCenter_PolarP.Hy);
            }


        }

        //花纹块正常左上轨迹
        public static void OfflineCalculateTrajectory_PatternUL(ref WorkPieceInfo workP)
        {
            double radius = Math.Sqrt(Pow(workP.CenterCir.X - workP.OfflineUpLeftPattern_OrthogonalP.X) + Pow(workP.CenterCir.Y - workP.OfflineUpLeftPattern_OrthogonalP.Y));
            workP.OfflinePartternPrFirst = radius;
            workP.OfflineULPartternCoordinates = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, workP.OfflineUpLeftPattern_OrthogonalP.X, workP.OfflineUpLeftPattern_OrthogonalP.Y, workP.OfflineSavePartternPrFirstStart.Y2, workP.OfflineUpLeftPattern_PolarP.Hy);

        }


        //花纹块正常右下

        public static void OfflineCalculateTrajectory_PatternDR(ref WorkPieceInfo workP)
        {
            double radius = Math.Sqrt(Pow(workP.CenterCir.X - workP.OfflineDownRightPattern_OrthogonalP.X) + Pow(workP.CenterCir.Y - workP.OfflineDownRightPattern_OrthogonalP.Y));
            workP.OfflinePartternPrSecond = radius;
            workP.OfflineDRPartternCoordinates = GetTrajectoryF(workP.CenterCir.X, workP.CenterCir.Y, workP.OfflineDownRightPattern_OrthogonalP.X, workP.OfflineDownRightPattern_OrthogonalP.Y, workP.OfflineSavePartternPrSecondStart.Y2, workP.OfflineDownRightPattern_PolarP.Hy);

        }


        //花纹块特殊轨迹

        public static void OfflineCalculateTrajectory_PatternSp(ref WorkPieceInfo workP)
        {

            List<Tuple<OrthogonalCoor, PolarCoor, double>> coordinates = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();

            double x1, y1, xb;

            if (workP.OfflineIsSmall)
            {
                x1 = workP.OfflineSPSmallPattern_OrthogonalP.X;
                y1 = workP.OfflineSPSmallPattern_OrthogonalP.Y;
                xb = workP.OfflineSPSmall_Pattern_PolarP.Hy;
            }
            else
            {
                x1 = workP.OfflineSPBigPattern_OrthogonalP.X;
                y1 = workP.OfflineSPBigPattern_OrthogonalP.Y;
                xb = workP.OfflineSPBig_Pattern_PolarP.Hy;
            }
            double radius = Math.Sqrt(Pow(workP.CenterCir.X - x1) + Pow(workP.CenterCir.Y - y1));
            workP.OfflinePartternPrThird = radius;
            workP.OfflineSPPartternCoordinates = GetTrajectory(workP.CenterCir.X, workP.CenterCir.Y, x1, y1, workP.OfflineSavePartternPrThirdStart.Y2, xb);

        }

        #endregion


        #region ptdata
        public static void RYPTData(List<Tuple<OrthogonalCoor, PolarCoor, double>> tuples, ref WorkPieceInfo workP, bool isCB, double radiusF)
        {
            workP.RptData.Clear();
            workP.YptData.Clear();
            workP.OptData.Clear();
            workP.RptDataF.Clear();
            workP.YptDataF.Clear();
            workP.OptDataF.Clear();

            Tuple<double, double, double> tuple = MotionTimeS(isCB, radiusF);
            workP.AllTime = tuple.Item3;
            workP.OVel = (360.0 + CommonMethod.flowPrmInfo.MarginAngle) * 1000 / workP.AllTime;
            double firstR = tuples[0].Item2.R;
            double firstY = tuples[0].Item2.Hy;
            double addR = 0, addY = 0, addtime = 0, addO = 0;
            double addRF = 0, addYF = 0, addtimeF = 0, addOF = 0;
            int dir = workP.OFollowDir ? 1 : -1;

            for (int i = 1; i < tuples.Count; i++)
            {
                addtime += tuple.Item2;
                addR += tuples[i].Item2.R - tuples[i - 1].Item2.R;
                addY += tuples[i].Item2.Hy - tuples[i - 1].Item2.Hy;
                addO += unit * dir;

                if (i == 1)
                {
                    workP.RptData.Add(new Tuple<double, double, short>(addR, addtime, PT_SEGMENT_NORMAL));
                    workP.YptData.Add(new Tuple<double, double, short>(addY, addtime, PT_SEGMENT_NORMAL));
                    workP.OptData.Add(new Tuple<double, double, short>(addO, addtime, PT_SEGMENT_NORMAL));
                }
                else if (i > 1 && i < tuples.Count - 1)
                {
                    workP.RptData.Add(new Tuple<double, double, short>(addR, addtime, PT_SEGMENT_EVEN));
                    workP.YptData.Add(new Tuple<double, double, short>(addY, addtime, PT_SEGMENT_EVEN));
                    workP.OptData.Add(new Tuple<double, double, short>(addO, addtime, PT_SEGMENT_EVEN));
                }
                else
                {
                    workP.RptData.Add(new Tuple<double, double, short>(addR, addtime, PT_SEGMENT_STOP));
                    workP.YptData.Add(new Tuple<double, double, short>(addY, addtime, PT_SEGMENT_STOP));
                    workP.OptData.Add(new Tuple<double, double, short>(addO, addtime, PT_SEGMENT_STOP));
                }
            }

            //workP.DownSideCoordinatesF = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();
            List<Tuple<OrthogonalCoor, PolarCoor, double>>  CoordinatesF = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();
            foreach (var item in tuples)
            {
                CoordinatesF.Add(item);
            }
            CoordinatesF.Reverse();

            for (int i = 1; i < CoordinatesF.Count; i++)
            {
                addtimeF += tuple.Item2;
                addRF += CoordinatesF[i].Item2.R - CoordinatesF[i - 1].Item2.R;
                addYF += CoordinatesF[i].Item2.Hy - CoordinatesF[i - 1].Item2.Hy;
                addOF += unit * (-dir);

                if (i == 1)
                {
                    workP.RptDataF.Add(new Tuple<double, double, short>(addRF, addtimeF, PT_SEGMENT_NORMAL));
                    workP.YptDataF.Add(new Tuple<double, double, short>(addYF, addtimeF, PT_SEGMENT_NORMAL));
                    workP.OptDataF.Add(new Tuple<double, double, short>(addOF, addtimeF, PT_SEGMENT_NORMAL));
                }
                else if (i > 1 && i < tuples.Count - 1)
                {
                    workP.RptDataF.Add(new Tuple<double, double, short>(addRF, addtimeF, PT_SEGMENT_EVEN));
                    workP.YptDataF.Add(new Tuple<double, double, short>(addYF, addtimeF, PT_SEGMENT_EVEN));
                    workP.OptDataF.Add(new Tuple<double, double, short>(addOF, addtimeF, PT_SEGMENT_EVEN));
                }
                else
                {
                    workP.RptDataF.Add(new Tuple<double, double, short>(addRF, addtimeF, PT_SEGMENT_STOP));
                    workP.YptDataF.Add(new Tuple<double, double, short>(addYF, addtimeF, PT_SEGMENT_STOP));
                    workP.OptDataF.Add(new Tuple<double, double, short>(addOF, addtimeF, PT_SEGMENT_STOP));
                }
            }



        }     
        #endregion


        #region time v
        public static Tuple<double, double> MotionTime(bool isCB, double radiusF)
        {
            double linearSpeed;
            if (isCB)
                linearSpeed = CommonMethod.flowPrmInfo.SideLineS;
            else
                linearSpeed = CommonMethod.flowPrmInfo.PatternLineS;

            double circumference = 2 * Math.PI * radiusF;

            double totalTime = circumference / linearSpeed * 1000;

            double segmentDistance = circumference / NumPoints;

            double timeStep = totalTime / NumPoints;

            double interpolationSpeed = segmentDistance / timeStep;

            return new Tuple<double, double>(totalTime, interpolationSpeed);
        }

        public static Tuple<double, double, double> MotionTimeS(bool isCB, double radiusF)
        {
            double linearSpeed;
            if (isCB)
            {
                if (CommonMethod.CBSuper)
                    linearSpeed = CommonMethod.flowPrmInfo.CBSuperLine;
                else
                    linearSpeed = CommonMethod.flowPrmInfo.SideLineS;
            }               
            else
            {
                if (CommonMethod.HWKSuper)
                    linearSpeed = CommonMethod.flowPrmInfo.HWKSuperLine;
                else
                    linearSpeed = CommonMethod.flowPrmInfo.PatternLineS;
            };

            double circumference = 2 * Math.PI * radiusF;
            circumference = circumference + circumference / 360.0 * CommonMethod.flowPrmInfo.MarginAngle;
            double totalTime = circumference / linearSpeed * 1000;

            double segmentDistance = circumference / NumPoints;

            double timeStep = totalTime / NumPoints;


            return new Tuple<double, double, double>(segmentDistance, timeStep, totalTime);
        }



        public static double VectorVelocityCal(bool isCB, double radius, short ul = 0)
        {
            double linearSpeed;
            if (isCB)
                linearSpeed = CommonMethod.flowPrmInfo.SideLineS;
            else
                linearSpeed = CommonMethod.flowPrmInfo.PatternLineS;

            //这里计算R1 Y2 O的矢量速度   线速度为mm/s单位
            double asp = (linearSpeed / radius) * 180 / Math.PI;

            double x = linearSpeed * CommonMethod.axisPrmInfo_R1.PosScale;
            double y = linearSpeed * CommonMethod.axisPrmInfo_Y2.PosScale;
            double z = asp * CommonMethod.axisPrmInfo_O.PosScale;

            double s = (CommonMethod.axisPrmInfo_R1.PosScale + CommonMethod.axisPrmInfo_Y2.PosScale + CommonMethod.axisPrmInfo_O.PosScale) / 3.0;
            //double v = linearSpeed / 10.0 * 128.6;
            double v;
            if (isCB)
            {
                v = (0.59 * linearSpeed / 10.0 * (128.6 * radius / 250)) * radius / 305;
            }
            else
            {
                v = linearSpeed / 10.0 * 60;
            }


            return v;
        }



        //跟随单位位移
        public static void FollowCal(ref WorkPieceInfo workP)
        {
            double angle = 360 + ExtraDegrees;

            workP.OFollowAngle = angle / NumPoints;

        }

        #endregion




        #region Common
        private static List<Tuple<OrthogonalCoor, PolarCoor, double>> GetTrajectory(double cX, double cY, double sX, double sY, double y2, double hy)
        {
            List<Tuple<OrthogonalCoor, PolarCoor, double>> coordinates = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();
            List<Tuple<Tuple<double, double>, Tuple<double, double>, double>> coors = CalculateTrajectory.Trajectory(cX, cY, sX, sY, ExtraDegrees, NumPoints, y2, hy);

            foreach (var item in coors)
            {
                coordinates.Add(new Tuple<OrthogonalCoor, PolarCoor, double>(new OrthogonalCoor(item.Item1.Item1, item.Item1.Item2), new PolarCoor(item.Item2.Item1, item.Item2.Item2), item.Item3));
            }

            return coordinates;

        }
        private static List<Tuple<OrthogonalCoor, PolarCoor, double>> GetTrajectoryF(double cX, double cY, double sX, double sY, double y2, double hy)
        {
            List<Tuple<OrthogonalCoor, PolarCoor, double>> coordinates = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();
            List<Tuple<Tuple<double, double>, Tuple<double, double>, double>> coors = CalculateTrajectory.TrajectoryF(cX, cY, sX, sY, ExtraDegrees, NumPoints, y2, hy);

            foreach (var item in coors)
            {
                coordinates.Add(new Tuple<OrthogonalCoor, PolarCoor, double>(new OrthogonalCoor(item.Item1.Item1, item.Item1.Item2), new PolarCoor(item.Item2.Item1, item.Item2.Item2), item.Item3));
            }

            return coordinates;

        }
        private static List<Tuple<OrthogonalCoor, PolarCoor, double>> GetSideTrajectory(double cX, double cY, double sX, double sY, double y2, double hy)
        {
            List<Tuple<OrthogonalCoor, PolarCoor, double>> coordinates = new List<Tuple<OrthogonalCoor, PolarCoor, double>>();
            List<Tuple<Tuple<double, double>, Tuple<double, double>, double>> coors = CalculateTrajectory.Trajectory(cX, cY, sX, sY, 0, 720, y2, hy);

            foreach (var item in coors)
            {
                coordinates.Add(new Tuple<OrthogonalCoor, PolarCoor, double>(new OrthogonalCoor(item.Item1.Item1, item.Item1.Item2), new PolarCoor(item.Item2.Item1, item.Item2.Item2), item.Item3));
            }

            return coordinates;

        }
        private static double Pow(double m)
        {
            return Math.Pow(m, 2);
        }
        private static double AR(double degrees)
        {
            return degrees * Math.PI / 180;
        }
        private static double RA(double radians)
        {
            return radians * 180 / Math.PI;
        }
        // 计算θ 从第四象限（0-90度）和第一象限（90-180度）方向
        private static double CalculateTheta(double x, double y)
        {
            if (x < 0)
                // X坐标为负
                return RA(Math.Atan2(y, Math.Abs(x))); // 转换到0-180度范围
            else
                // X坐标为正
                return 90 + (90 - RA(Math.Atan2(y, x)));
        }


        private static Tuple<double, double> DownCB_R3Rota_OffYZ(double rotAngle)
        {
            double shortb = N;
            double longb = DisFocus;
            double hy = Math.Sqrt(Pow(shortb) + Pow(longb));

            double y, z;

            double includedAngle = 2 * RA(Math.Atan(shortb / longb));
            if (rotAngle == 0) { y = 0; z = 0; }
            else if (rotAngle > 0 && rotAngle < includedAngle)
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = x1 - (90 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = -Math.Cos(AR(a1)) * b1;
                z = Math.Sin(AR(a1)) * b1;
            }
            else if (rotAngle > includedAngle)
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = 90 - Math.Atan(shortb / longb) * 180 / Math.PI - x1;
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = -Math.Cos(AR(a1)) * b1;
                z = -Math.Sin(AR(a1)) * b1;
            }
            else if (rotAngle < 0)
            {
                double x1 = (180 - Math.Abs(rotAngle)) / 2;
                double a1 = 90 - (x1 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(Math.Abs(rotAngle) / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = -Math.Sin(AR(a1)) * b1;
            }
            else
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = x1 - (90 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = 0;
            }

            return new Tuple<double, double>(y, z);
        }

        private static Tuple<double, double> UpCB_R3Rota_OffYZ(double rotAngle)
        {
            double shortb = N;
            double longb = DisFocus;
            double hy = Math.Sqrt(Pow(shortb) + Pow(longb));

            double y, z;

            double includedAngle = 2 * RA(Math.Atan(shortb / longb));
            if (rotAngle == 0) { y = 0; z = 0; }
            else if (rotAngle > 0 && rotAngle < includedAngle)
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = x1 - (90 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = -Math.Sin(AR(a1)) * b1;
            }
            else if (rotAngle > includedAngle)
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = 90 - Math.Atan(shortb / longb) * 180 / Math.PI - x1;
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = Math.Sin(AR(a1)) * b1;
            }
            else if (rotAngle < 0)
            {
                double x1 = (180 - Math.Abs(rotAngle)) / 2;
                double a1 = 90 - (x1 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(Math.Abs(rotAngle) / 2));
                y = -Math.Cos(AR(a1)) * b1;
                z = Math.Sin(AR(a1)) * b1;
            }
            else
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = x1 - (90 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = 0;
            }

            return new Tuple<double, double>(y, z);
        }

        private static Tuple<double, double> Pattern_R3Rota_OffYZ(double rotAngle, short typeN)
        {
            double shortb = N;
            double longb = 0;

            switch (typeN)
            {
                case 0:
                    longb = DisFocus;
                    break;
                case 1:
                    longb = DisFocus + CommonMethod.flowPrmInfo.ULZoomOff;
                    break;
                case 2:
                    longb = DisFocus + CommonMethod.flowPrmInfo.DRZoomOff;
                    break;
                case 3:
                    longb = DisFocus + CommonMethod.flowPrmInfo.OfflineULZoomOff;
                    break;
                case 4:
                    longb = DisFocus + CommonMethod.flowPrmInfo.OfflineDRZoomOff;
                    break;
                default:
                    longb = DisFocus;
                    break;
            }
            double hy = Math.Sqrt(Pow(shortb) + Pow(longb));

            double y, z;

            double includedAngle = 2 * RA(Math.Atan(shortb / longb));
            if (rotAngle == 0) { y = 0; z = 0; }
            else if (rotAngle > 0 && rotAngle < includedAngle)
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = x1 - (90 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = Math.Sin(AR(a1)) * b1;
            }
            else if (rotAngle > includedAngle)
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = 90 - Math.Atan(shortb / longb) * 180 / Math.PI - x1;
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = -Math.Sin(AR(a1)) * b1;
            }
            else if (rotAngle < 0)
            {
                double x1 = (180 - Math.Abs(rotAngle)) / 2;
                double a1 = 90 - (x1 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(Math.Abs(rotAngle) / 2));
                y = -Math.Cos(AR(a1)) * b1;
                z = -Math.Sin(AR(a1)) * b1;
            }
            else
            {
                double x1 = (180 - rotAngle) / 2;
                double a1 = x1 - (90 - Math.Atan(shortb / longb) * 180 / Math.PI);
                double b1 = 2 * hy * Math.Sin(AR(rotAngle / 2));
                y = Math.Cos(AR(a1)) * b1;
                z = 0;
            }

            return new Tuple<double, double>(y, z);
        }


        //求上下点n边映射的同心圆的n边长度
        private static double LR_CalculateN(bool isLeft)
        {
            //左朝上  右朝下(朝下的角度有三种算法)
            double n;
            double xb = Math.Abs(Math.Sqrt(Pow(M) + Pow(DisFocus)));
            double angle1 = RA(Math.Asin(M / xb));
            if (isLeft)
            {
                double angle2 = angle1 + CommonMethod.workPieceInfo.UA;
                n = xb * Math.Cos(AR(angle2));

            }
            else
            {
                if (CommonMethod.workPieceInfo.DA < angle1)
                {
                    double angle2 = angle1 - CommonMethod.workPieceInfo.DA;
                    n = xb * Math.Cos(AR(angle2));
                }
                else if (CommonMethod.workPieceInfo.DA == angle1)
                {
                    n = xb;
                }
                else
                {
                    double angle2 = 90 - CommonMethod.workPieceInfo.DA + angle1;
                    n = xb * Math.Sin(AR(angle2));
                }

            }

            return n;
        }

        private static double LR_CalculateN2(bool isLeft, short typeN)
        {
            double dis = 0;
            switch (typeN)
            {
                case 0:
                    dis = DisFocus;
                    break;
                case 1:
                    dis = DisFocus + CommonMethod.flowPrmInfo.ULZoomOff;
                    break;
                case 2:
                    dis = DisFocus + CommonMethod.flowPrmInfo.DRZoomOff;
                    break;
                case 3:
                    dis = DisFocus + CommonMethod.flowPrmInfo.OfflineULZoomOff;
                    break;
                case 4:
                    dis = DisFocus + CommonMethod.flowPrmInfo.OfflineDRZoomOff;
                    break;
                default:
                    dis = DisFocus;
                    break;
            }
            //左朝上  右朝下
            double newn;
            double xxb = Math.Abs(Math.Sqrt(Pow(N) + Pow(dis)));
            double xxbAn = RA(Math.Asin(dis / xxb));
            double xxbAd = RA(Math.Asin(N / xxb));

            if (isLeft)
            {
                double angle2 = 90 - xxbAn + CommonMethod.workPieceInfo.UA;
                newn = xxb * Math.Cos(AR(angle2));

            }
            else
            {
                if (CommonMethod.workPieceInfo.DA < xxbAd)
                {
                    double angle2 = xxbAd - CommonMethod.workPieceInfo.DA;
                    newn = xxb * Math.Cos(AR(angle2));
                }
                else if (CommonMethod.workPieceInfo.DA == xxbAd)
                {
                    newn = xxb;
                }
                else
                {
                    double angle2 = CommonMethod.workPieceInfo.DA - xxbAd;
                    newn = xxb * Math.Cos(AR(angle2));
                }

            }

            return newn;
        }

        //求n外端点角度
        private static double L_CalculateNR_A(double n, WorkPieceInfo workP)
        {
            double angle1 = RA(Math.Atan(n / M));
            double angle2 = 90 - angle1;
            double xb = Math.Sqrt(Pow(n) + Pow(M));
            double by = workP.UpPattern_PolarP.Hy;
            double c = Math.Abs(Math.Sqrt(Pow(xb) + Pow(by) - 2 * xb * by * Math.Cos(AR(angle2))));

            double angle = RA(Math.Acos((Pow(by) + Pow(c) - Pow(xb)) / (2 * by * c)));

            return angle + workP.UpPattern_PolarP.R;
        }

        private static double R_CalculateNR_A(double n, WorkPieceInfo workP)
        {
            double angle1 = RA(Math.Atan(n / M));
            double angle2 = 90 - angle1;
            double xb = Math.Sqrt(Pow(n) + Pow(M));
            double by = workP.DownPattern_PolarP.Hy;
            double c = Math.Abs(Math.Sqrt(Pow(xb) + Pow(by) - 2 * xb * by * Math.Cos(AR(angle2))));

            double angle = RA(Math.Acos((Pow(by) + Pow(c) - Pow(xb)) / (2 * by * c)));

            return angle + workP.DownPattern_PolarP.R;
        }

        //求n外端点角度
        private static double OfflineL_CalculateNR_A(double n, WorkPieceInfo workP)
        {
            double angle1 = RA(Math.Atan(n / M));
            double angle2 = 90 - angle1;
            double xb = Math.Sqrt(Pow(n) + Pow(M));
            double by = workP.OfflineUpPattern_PolarP.Hy;
            double c = Math.Abs(Math.Sqrt(Pow(xb) + Pow(by) - 2 * xb * by * Math.Cos(AR(angle2))));

            double angle = RA(Math.Acos((Pow(by) + Pow(c) - Pow(xb)) / (2 * by * c)));

            return angle + workP.OfflineUpPattern_PolarP.R;
        }

        private static double OfflineR_CalculateNR_A(double n, WorkPieceInfo workP)
        {
            double angle1 = RA(Math.Atan(n / M));
            double angle2 = 90 - angle1;
            double xb = Math.Sqrt(Pow(n) + Pow(M));
            double by = workP.OfflineDownPattern_PolarP.Hy;
            double c = Math.Abs(Math.Sqrt(Pow(xb) + Pow(by) - 2 * xb * by * Math.Cos(AR(angle2))));

            double angle = RA(Math.Acos((Pow(by) + Pow(c) - Pow(xb)) / (2 * by * c)));

            return angle + workP.OfflineDownPattern_PolarP.R;
        }

        private static Tuple<double, double> Polar_Orthogonal(PolarCoor polarCoor)
        {
            double x, y;
            if (polarCoor.R < 90)
            {
                x = -polarCoor.Hy * Math.Cos(AR(polarCoor.R));
                y = polarCoor.Hy * Math.Sin(AR(polarCoor.R));
            }
            else
            {
                x = polarCoor.Hy * Math.Cos(AR(180 - polarCoor.R));
                y = polarCoor.Hy * Math.Sin(AR(180 - polarCoor.R));
            }

            return new Tuple<double, double>(x, y);
        }

        //左夹角实现
        private static void LAngleRealize(ref WorkPieceInfo workP)
        {
            // 圆的参数  这里采用2方法，计算新n
            double cx = workP.CenterCir.X;
            double cy = workP.CenterCir.Y;
            double r = workP.TreadDiameter / 2.0;

            double m = M;
            double n = LR_CalculateN2(true, 1);

            double initialMx = workP.UpPattern_OrthogonalP.X;
            double initialMy = workP.UpPattern_OrthogonalP.Y;

            double angle = L_CalculateNR_A(n, workP);
            double radian = (360 - angle) * Math.PI / 180;

            double x1 = cx + r * Math.Cos(AR(workP.UpPattern_PolarP.R));//cx + r * Math.Cos(radian)
            double y1 = cy - r * Math.Sin(AR(workP.UpPattern_PolarP.R));// cy + r * Math.Sin(radian)

            double deltaX = initialMx - x1;
            double deltaY = initialMy - y1;

            double cosTheta = Math.Cos(AR(workP.LA));
            double sinTheta = Math.Sin(AR(workP.LA));

            double deltaX_new = cosTheta * deltaX - sinTheta * deltaY;
            double deltaY_new = sinTheta * deltaX + cosTheta * deltaY;

            double newMx = x1 + deltaX_new;
            double newMy = y1 + deltaY_new;
            //ab夹角
            double c = Math.Sqrt(Pow(x1) + Pow(y1));
            double a = Math.Sqrt(Pow(m) + Pow(n));
            double b = Math.Sqrt(Pow(newMx) + Pow(newMy));

            double abAngle = RA(Math.Acos((Pow(a) + Pow(b) - Pow(c)) / (2 * a * b)));
            double sangle = RA(Math.Acos(m / a));

            double newAngle = abAngle + sangle - 90;

            double ra;
            if (newMx < 0)
                ra = Math.Abs(RA(Math.Atan(newMy / newMx)));
            else
                ra = 180 - RA(Math.Atan(newMy / newMx));


            //计算当前pr
            double pr = Math.Sqrt(Pow(cx - newMx) + Pow(cy - newMy));

            workP.UpLeftPattern_PolarP.R = ra;
            workP.UpLeftPattern_PolarP.Hy = b;
            workP.UpLeftPattern_OrthogonalP.X = newMx;
            workP.UpLeftPattern_OrthogonalP.Y = newMy;
            workP.UpLeftPattern_Pr = pr;
            workP.LeftPattern_OAngle = -Math.Abs(newAngle);



        }


        //右夹角实现
        private static void RAngleRealize(ref WorkPieceInfo workP)
        {
            double cx = workP.CenterCir.X;
            double cy = workP.CenterCir.Y;
            double r = workP.TreadDiameter / 2.0;

            double m = M;
            double n = LR_CalculateN2(false, 2);

            double initialMx = workP.DownPattern_OrthogonalP.X;
            double initialMy = workP.DownPattern_OrthogonalP.Y;

            double angle = R_CalculateNR_A(n, workP);
            double radian = (360 - angle) * Math.PI / 180;

            double x1 = cx + r * Math.Cos(AR(workP.DownPattern_PolarP.R));//cx + r * Math.Cos(radian)
            double y1 = cy - r * Math.Sin(AR(workP.DownPattern_PolarP.R));// cy + r * Math.Sin(radian)

            double deltaX = initialMx - x1;
            double deltaY = initialMy - y1;

            double cosTheta = Math.Cos(AR(-workP.RA));
            double sinTheta = Math.Sin(AR(-workP.RA));

            double deltaX_new = cosTheta * deltaX - sinTheta * deltaY;
            double deltaY_new = sinTheta * deltaX + cosTheta * deltaY;

            double newMx = x1 + deltaX_new;
            double newMy = y1 + deltaY_new;

            //ab夹角
            double c = Math.Sqrt(Pow(x1) + Pow(y1));
            double a = Math.Sqrt(Pow(m) + Pow(n));
            double b = Math.Sqrt(Pow(newMx) + Pow(newMy));

            double abAngle = RA(Math.Acos((Pow(a) + Pow(b) - Pow(c)) / (2 * a * b)));
            double sangle = RA(Math.Acos(m / a));


            double ra;
            if (newMx < 0)
                ra = Math.Abs(RA(Math.Atan(newMy / newMx)));
            else
                ra = 180 - RA(Math.Atan(newMy / newMx));

            double newAngle;

            if (ra > workP.DownCenter_PolarP.R + workP.Start_Rmn_A)
            {
                newAngle = 90 - (sangle - abAngle);
            }
            else if (ra == workP.DownCenter_PolarP.R + workP.Start_Rmn_A)
            {
                newAngle = sangle;
            }
            else
            {
                newAngle = 90 - (sangle + abAngle);
            }

            double pr = Math.Sqrt(Pow(cx - newMx) + Pow(cy - newMy));

            workP.DownRightPattern_PolarP.R = ra;
            workP.DownRightPattern_PolarP.Hy = b;
            workP.DownRightPattern_OrthogonalP.X = newMx;
            workP.DownRightPattern_OrthogonalP.Y = newMy;
            workP.DownRightPattern_Pr = pr;
            workP.RightPattern_OAngle = Math.Abs(newAngle);



        }

        //OA计算  角度是按二维的正角度旋转寻找 方向为逆时针， 输出为正值，实际需要负转 这里已进行取反处理
        private static ComResult OAngleCalculate(ref WorkPieceInfo workP)
        {
            double xA = workP.SPBigPattern_OrthogonalP.X, yA = workP.SPBigPattern_OrthogonalP.Y;
            double xB = workP.SPN_OrthogonalP.X, yB = workP.SPN_OrthogonalP.Y;
            double xC = workP.CenterCir.X, yC = workP.CenterCir.Y;
            double R = workP.TreadDiameter / 2.0;

            double dBC = Math.Sqrt(Math.Pow(xB - xC, 2) + Math.Pow(yB - yC, 2));

            if (dBC > R)
                return ComResult.CreateFailResult("无法回到圆上");

            double relX = xB - xA;
            double relY = yB - yA;

            bool foundSolution = false;
            double thetaRad = 0;

            for (int i = 0; i < 36000; i++)
            {
                thetaRad = Math.PI * i / 18000.0;

                double rotatedX = xA + relX * Math.Cos(thetaRad) - relY * Math.Sin(thetaRad);
                double rotatedY = yA + relX * Math.Sin(thetaRad) + relY * Math.Cos(thetaRad);

                double distanceToCenter = Math.Sqrt(Math.Pow(rotatedX - xC, 2) + Math.Pow(rotatedY - yC, 2));

                if (Math.Abs(distanceToCenter - R) < 3)//3mm
                {
                    foundSolution = true;
                    break;
                }
            }

            if (foundSolution)
            {
                double thetaDeg = thetaRad * 180 / Math.PI;
                workP.BigPattern_OAngle = -thetaDeg;
            }
            else
            {
                return ComResult.CreateFailResult("无法回到圆上");
            }

            return ComResult.CreateSuccessResult("可回到圆上，已计算完成");

        }

        public static void LCenterCal(int type=1)
        {
            double by = DisRY + CommonMethod.CentringCoordinate.Y1 + CommonMethod.CentringCoordinate.Y2;


            double x, y;
            if (CommonMethod.CentringCoordinate.R <= 90)
            {
                x = -by * Math.Cos(AR(CommonMethod.CentringCoordinate.R));
                y = by * Math.Sin(AR(CommonMethod.CentringCoordinate.R));
            }
            else
            {
                x = by * Math.Cos(AR(180 - CommonMethod.CentringCoordinate.R));
                y = by * Math.Sin(AR(180 - CommonMethod.CentringCoordinate.R));
            }
           
            double newx = 0, newy = 0;

            switch (type)
            {
                case 1:
                    newx = CommonMethod.CentringCoordinate.R <= 90 ? x + CommonMethod.flowPrmInfo.CenterCoorXYOff_C.X : x + CommonMethod.flowPrmInfo.CenterCoorXYOff_CR.X;
                    newy = CommonMethod.CentringCoordinate.R <= 90 ? y + CommonMethod.flowPrmInfo.CenterCoorXYOff_C.Y : y + CommonMethod.flowPrmInfo.CenterCoorXYOff_CR.Y;
                    break;
                case 2://H1
                    newx = CommonMethod.CentringCoordinate.R <= 90 ? x + CommonMethod.flowPrmInfo.CenterCoorXYOff_LH1.X : x + CommonMethod.flowPrmInfo.CenterCoorXYOff_RH1.X;
                    newy = CommonMethod.CentringCoordinate.R <= 90 ? y + CommonMethod.flowPrmInfo.CenterCoorXYOff_LH1.Y : y + CommonMethod.flowPrmInfo.CenterCoorXYOff_RH1.Y;
                    break;
                case 3://H2
                    newx = CommonMethod.CentringCoordinate.R <= 90 ? x + CommonMethod.flowPrmInfo.CenterCoorXYOff_LH2.X : x + CommonMethod.flowPrmInfo.CenterCoorXYOff_RH2.X;
                    newy = CommonMethod.CentringCoordinate.R <= 90 ? y + CommonMethod.flowPrmInfo.CenterCoorXYOff_LH2.Y : y + CommonMethod.flowPrmInfo.CenterCoorXYOff_RH2.Y;
                    break;                
                default:
                    newx = 0; newy = 0;
                    break;
            }


            double newby = Math.Abs(Math.Sqrt(Pow(newx) + Pow(newy)));
            double newr = newx <= 0 ? RA(Math.Atan2(newy, -newx)) : 180 - RA(Math.Atan2(newy, newx));
            double newy2 = CommonMethod.CentringCoordinate.Y2 + (newby - by);

            double hby = newby + CommonMethod.flowPrmInfo.H_O_HorDist;
            double xb = Math.Sqrt(Pow(N) + Pow(M));
            double a1 = 90 - RA(Math.Asin(N / xb));

            double c = xb * Math.Cos(AR(a1)) + Math.Abs(Math.Sqrt(Pow(hby) - Pow(xb) * Pow(Math.Sin(AR(a1)))));

            double ang = RA(Math.Acos((Pow(c) + Pow(hby) - Pow(xb)) / (2 * c * hby)));

            CommonMethod.CentringCoordinate_L.R = newr - ang;
            CommonMethod.CentringCoordinate_L.Y2 = newy2 + (c - newby);
            CommonMethod.CentringCoordinate_L.Y1 = CommonMethod.CentringCoordinate.Y1;
            CommonMethod.CentringCoordinate_L.Z = CommonMethod.CentringCoordinate.Z;


        }
        #endregion

        #region offline common
        //左夹角实现
        private static void OfflineLAngleRealize(ref WorkPieceInfo workP)
        {
            double cx = workP.CenterCir.X;
            double cy = workP.CenterCir.Y;
            double r = workP.TreadDiameter / 2.0;

            double m = M;
            double n = LR_CalculateN2(true, 3);

            double initialMx = workP.OfflineUpPattern_OrthogonalP.X;
            double initialMy = workP.OfflineUpPattern_OrthogonalP.Y;

            double angle = OfflineL_CalculateNR_A(n, workP);
            double radian = (360 - angle) * Math.PI / 180;

            double x1 = cx + r * Math.Cos(AR(workP.OfflineUpPattern_PolarP.R));//cx + r * Math.Cos(radian)
            double y1 = cy - r * Math.Sin(AR(workP.OfflineUpPattern_PolarP.R));// cy + r * Math.Sin(radian)
            double deltaX = initialMx - x1;
            double deltaY = initialMy - y1;
            double cosTheta = Math.Cos(AR(workP.LA));
            double sinTheta = Math.Sin(AR(workP.LA));

            double deltaX_new = cosTheta * deltaX - sinTheta * deltaY;
            double deltaY_new = sinTheta * deltaX + cosTheta * deltaY;

            double newMx = x1 + deltaX_new;
            double newMy = y1 + deltaY_new;

            //ab夹角
            double c = Math.Sqrt(Pow(x1) + Pow(y1));
            double a = Math.Sqrt(Pow(m) + Pow(n));
            double b = Math.Sqrt(Pow(newMx) + Pow(newMy));

            double abAngle = RA(Math.Acos((Pow(a) + Pow(b) - Pow(c)) / (2 * a * b)));
            double sangle = RA(Math.Acos(m / a));

            double newAngle = abAngle + sangle - 90;

            double ra;
            if (newMx < 0)
                ra = Math.Abs(RA(Math.Atan(newMy / newMx)));
            else
                ra = 180 - RA(Math.Atan(newMy / newMx));


            double pr = Math.Sqrt(Pow(cx - newMx) + Pow(cy - newMy));

            workP.OfflineUpLeftPattern_PolarP.R = ra;
            workP.OfflineUpLeftPattern_PolarP.Hy = b;
            workP.OfflineUpLeftPattern_OrthogonalP.X = newMx;
            workP.OfflineUpLeftPattern_OrthogonalP.Y = newMy;
            workP.OfflineUpLeftPattern_Pr = pr;
            workP.OfflineLeftPattern_OAngle = -Math.Abs(newAngle);



        }


        //右夹角实现
        private static void OfflineRAngleRealize(ref WorkPieceInfo workP)
        {
            double cx = workP.CenterCir.X;
            double cy = workP.CenterCir.Y;
            double r = workP.TreadDiameter / 2.0;

            double m = M;
            double n = LR_CalculateN2(false, 4);
            double initialMx = workP.OfflineDownPattern_OrthogonalP.X;
            double initialMy = workP.OfflineDownPattern_OrthogonalP.Y;
            double angle = OfflineR_CalculateNR_A(n, workP);
            double radian = (360 - angle) * Math.PI / 180;

            double x1 = cx + r * Math.Cos(AR(workP.OfflineDownPattern_PolarP.R));//cx + r * Math.Cos(radian)
            double y1 = cy - r * Math.Sin(AR(workP.OfflineDownPattern_PolarP.R));// cy + r * Math.Sin(radian)

            double deltaX = initialMx - x1;
            double deltaY = initialMy - y1;

            double cosTheta = Math.Cos(AR(-workP.RA));
            double sinTheta = Math.Sin(AR(-workP.RA));

            double deltaX_new = cosTheta * deltaX - sinTheta * deltaY;
            double deltaY_new = sinTheta * deltaX + cosTheta * deltaY;

            double newMx = x1 + deltaX_new;
            double newMy = y1 + deltaY_new;

            //ab夹角
            double c = Math.Sqrt(Pow(x1) + Pow(y1));
            double a = Math.Sqrt(Pow(m) + Pow(n));
            double b = Math.Sqrt(Pow(newMx) + Pow(newMy));

            double abAngle = RA(Math.Acos((Pow(a) + Pow(b) - Pow(c)) / (2 * a * b)));
            double sangle = RA(Math.Acos(m / a));


            double ra;
            if (newMx < 0)
                ra = Math.Abs(RA(Math.Atan(newMy / newMx)));
            else
                ra = 180 - RA(Math.Atan(newMy / newMx));

            double newAngle;

            if (ra > workP.DownCenter_PolarP.R + workP.Start_Rmn_A)
            {
                newAngle = 90 - (sangle - abAngle);
            }
            else if (ra == workP.DownCenter_PolarP.R + workP.Start_Rmn_A)
            {
                newAngle = sangle;
            }
            else
            {
                newAngle = 90 - (sangle + abAngle);
            }

            double pr = Math.Sqrt(Pow(cx - newMx) + Pow(cy - newMy));

            workP.OfflineDownRightPattern_PolarP.R = ra;
            workP.OfflineDownRightPattern_PolarP.Hy = b;
            workP.OfflineDownRightPattern_OrthogonalP.X = newMx;
            workP.OfflineDownRightPattern_OrthogonalP.Y = newMy;
            workP.OfflineDownRightPattern_Pr = pr;
            workP.OfflineRightPattern_OAngle = Math.Abs(newAngle);



        }

        //OA计算  角度是按二维的正角度旋转寻找 方向为逆时针， 输出为正值，实际需要负转 这里已进行取反处理
        private static ComResult OfflineOAngleCalculate(ref WorkPieceInfo workP)
        {
            double xA = workP.OfflineSPBigPattern_OrthogonalP.X, yA = workP.OfflineSPBigPattern_OrthogonalP.Y;
            double xB = workP.OfflineSPN_OrthogonalP.X, yB = workP.OfflineSPN_OrthogonalP.Y;
            double xC = workP.CenterCir.X, yC = workP.CenterCir.Y;
            double R = workP.TreadDiameter / 2.0;

            double dBC = Math.Sqrt(Math.Pow(xB - xC, 2) + Math.Pow(yB - yC, 2));

            if (dBC > R)
                return ComResult.CreateFailResult("无法回到圆上");

            double relX = xB - xA;
            double relY = yB - yA;

            // 声明旋转角度变量
            bool foundSolution = false;
            double thetaRad = 0;

            for (int i = 0; i < 36000; i++)
            {
                thetaRad = Math.PI * i / 18000.0;

                double rotatedX = xA + relX * Math.Cos(thetaRad) - relY * Math.Sin(thetaRad);
                double rotatedY = yA + relX * Math.Sin(thetaRad) + relY * Math.Cos(thetaRad);

                double distanceToCenter = Math.Sqrt(Math.Pow(rotatedX - xC, 2) + Math.Pow(rotatedY - yC, 2));

                if (Math.Abs(distanceToCenter - R) < 3)//3mm
                {
                    foundSolution = true;
                    break;
                }
            }

            if (foundSolution)
            {
                double thetaDeg = thetaRad * 180 / Math.PI;
                workP.OfflineBigPattern_OAngle = -thetaDeg;
            }
            else
            {
                return ComResult.CreateFailResult("无法回到圆上");
            }

            return ComResult.CreateSuccessResult("可回到圆上，已计算完成");

        }
        #endregion


    }
}

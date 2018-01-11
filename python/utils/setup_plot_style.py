"""
Setup plotting style according to
https://twiki.cern.ch/twiki/bin/view/CMS/Internal/FigGuidelines
"""
import ROOT as r

def set_TDR_style():
    tdr_style = r.TStyle('tdr_style', 'Style for P-TDR')

    # canvas settings
    tdr_style.SetCanvasBorderMode(0)
    tdr_style.SetCanvasColor(r.kWhite)
    tdr_style.SetCanvasDefH(600) # height
    tdr_style.SetCanvasDefW(600) # width
    tdr_style.SetCanvasDefX(0) # position on screen
    tdr_style.SetCanvasDefY(0)

    tdr_style.SetPadBorderMode(0)
    # tdr_style.SetPadBorderSize(Width_t size = 1)
    tdr_style.SetPadColor(r.kWhite)
    tdr_style.SetPadGridX(False)
    tdr_style.SetPadGridY(False)
    tdr_style.SetGridColor(0)
    tdr_style.SetGridStyle(3)
    tdr_style.SetGridWidth(1)

    #For the frame:
    tdr_style.SetFrameBorderMode(0)
    tdr_style.SetFrameBorderSize(1)
    tdr_style.SetFrameFillColor(0)
    tdr_style.SetFrameFillStyle(0)
    tdr_style.SetFrameLineColor(1)
    tdr_style.SetFrameLineStyle(1)
    tdr_style.SetFrameLineWidth(1)

    # #For the histo:
    # #tdr_style.SetHistFillColor(1)
    # #tdr_style.SetHistFillStyle(0)
    # tdr_style.SetHistLineColor(1)
    # tdr_style.SetHistLineStyle(0)
    # tdr_style.SetHistLineWidth(1)
    # #tdr_style.SetLegoInnerR(Float_t rad = 0.5)
    # #tdr_style.SetNumberContours(Int_t number = 20)

    # tdr_style.SetEndErrorSize(2)
    # #tdr_style.SetErrorMarker(20)
    # #tdr_style.SetErrorX(0.)

    # tdr_style.SetMarkerStyle(20)

    # #For the fit/function:
    # tdr_style.SetOptFit(1)
    # tdr_style.SetFitFormat("5.4g")
    # tdr_style.SetFuncColor(2)
    # tdr_style.SetFuncStyle(1)
    # tdr_style.SetFuncWidth(1)

    #For the date:
    tdr_style.SetOptDate(0)
    # tdr_style.SetDateX(Float_t x = 0.01)
    # tdr_style.SetDateY(Float_t y = 0.01)

    # For the statistics box:
    tdr_style.SetOptFile(0)
    tdr_style.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
    tdr_style.SetStatColor(r.kWhite)
    tdr_style.SetStatFont(42)
    tdr_style.SetStatFontSize(0.025)
    tdr_style.SetStatTextColor(1)
    tdr_style.SetStatFormat("6.4g")
    tdr_style.SetStatBorderSize(1)
    tdr_style.SetStatH(0.1)
    tdr_style.SetStatW(0.15)
    # tdr_style.SetStatStyle(Style_t style = 1001)
    # tdr_style.SetStatX(Float_t x = 0)
    # tdr_style.SetStatY(Float_t y = 0)

    # Margins:
    tdr_style.SetPadTopMargin(0.05)
    tdr_style.SetPadBottomMargin(0.13)
    tdr_style.SetPadLeftMargin(0.16)
    tdr_style.SetPadRightMargin(0.02)

    # For the Global title:

    tdr_style.SetOptTitle(0)
    tdr_style.SetTitleFont(42)
    tdr_style.SetTitleColor(1)
    tdr_style.SetTitleTextColor(1)
    tdr_style.SetTitleFillColor(10)
    tdr_style.SetTitleFontSize(0.05)
    # tdr_style.SetTitleH(0) # Set the height of the title box
    # tdr_style.SetTitleW(0) # Set the width of the title box
    # tdr_style.SetTitleX(0) # Set the position of the title box
    # tdr_style.SetTitleY(0.985) # Set the position of the title box
    # tdr_style.SetTitleStyle(Style_t style = 1001)
    # tdr_style.SetTitleBorderSize(2)

    # For the axis titles:

    tdr_style.SetTitleColor(1, "XYZ")
    tdr_style.SetTitleFont(42, "XYZ")
    tdr_style.SetTitleSize(0.06, "XYZ")
    # tdr_style.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
    # tdr_style.SetTitleYSize(Float_t size = 0.02)
    tdr_style.SetTitleXOffset(0.9)
    tdr_style.SetTitleYOffset(1.25)
    # tdr_style.SetTitleOffset(1.1, "Y") # Another way to set the Offset

    # For the axis labels:

    tdr_style.SetLabelColor(1, "XYZ")
    tdr_style.SetLabelFont(42, "XYZ")
    tdr_style.SetLabelOffset(0.007, "XYZ")
    tdr_style.SetLabelSize(0.05, "XYZ")

    # For the axis:

    tdr_style.SetAxisColor(1, "XYZ")
    tdr_style.SetStripDecimals(True)
    tdr_style.SetTickLength(0.03, "XYZ")
    tdr_style.SetNdivisions(510, "XYZ")
    tdr_style.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
    tdr_style.SetPadTickY(1)

    # Change for log plots:
    tdr_style.SetOptLogx(0)
    tdr_style.SetOptLogy(0)
    tdr_style.SetOptLogz(0)

    # Postscript options:
    tdr_style.SetPaperSize(20. ,20.)
    # tdr_style.SetLineScalePS(Float_t scale = 3)
    # tdr_style.SetLineStyleString(Int_t i, const char* text)
    # tdr_style.SetHeaderPS(const char* header)
    # tdr_style.SetTitlePS(const char* pstitle)

    # tdr_style.SetBarOffset(Float_t baroff = 0.5)
    # tdr_style.SetBarWidth(Float_t barwidth = 0.5)
    # tdr_style.SetPaintTextFormat(const char* format = "g")
    # tdr_style.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
    # tdr_style.SetTimeOffset(Double_t toffset)
    # tdr_style.SetHistMinimumZero(true)

    tdr_style.SetHatchesLineWidth(5)
    tdr_style.SetHatchesSpacing(0.05)

    tdr_style.cd()


def get_pad_margins(pad):
    """Get all the necessary info from the pad and put it into a dict"""
    pad_info = {}
    pad_info['H'] = pad.GetWh()
    pad_info['W'] = pad.GetWw()
    pad_info['l'] = pad.GetLeftMargin()
    pad_info['r'] = pad.GetRightMargin()
    pad_info['t'] = pad.GetTopMargin()
    pad_info['b'] = pad.GetBottomMargin()
    pad_info['e'] = 0.25

    return pad_info


def add_auxiliary_info(pad, years, pos='right'):
    """Add the auxiliary information to the passed pad"""
    LUMINOSITY = {'2012': '19.7 fb^{-1} (8 TeV)',
                  '2016': '4.6 fb^{-1} (13 TeV)',
                  '2017': '42.4 fb^{-1} (13 TeV)'}

    # required text setup
    CMS_TEXT = 'CMS'
    CMS_TEXT_FONT = 61
    CMS_TEXT_SIZE = 0.75
    CMS_TEXT_OFFSET = 0.1

    EXTRA_TEXT = 'Preliminary'
    EXTRA_TEXT_FONT = 52
    EXTRA_TEXT_SIZE = 0.76 * CMS_TEXT_SIZE

    LUMI_TEXT_SIZE = 0.6
    LUMI_TEXT_OFFSET = 0.2

    REL_POS_X = 0.045
    REL_POS_Y = 0.065
    REL_EXTRA_DY = 1.1


    PAD = get_pad_margins(pad)

    pad.cd()
    # latex used to draw everything
    latex = r.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(r.kBlack)

    # lumi info
    lumi_text = ' + '.join([LUMINOSITY[y] for y in years])
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(LUMI_TEXT_SIZE * PAD['t'])
    latex.DrawLatex(1-PAD['r'], 1-PAD['t']+LUMI_TEXT_OFFSET*PAD['t'], lumi_text)

    if pos == 'right':
        pos_x = 1 - PAD['r'] - REL_POS_X * (1 - PAD['l'] - PAD['r'])
    if pos == 'left':
        pos_x = PAD['l'] + REL_POS_X * (1 - PAD['l'] - PAD['r'])

    pos_y = 1 - PAD['t'] - REL_POS_Y * (1 - PAD['t'] - PAD['b'])

    text_align = 11
    if pos == 'right':
        text_align = 31

    # cms text
    latex.SetTextFont(CMS_TEXT_FONT)
    latex.SetTextSize(CMS_TEXT_SIZE * PAD['t'])
    latex.SetTextAlign(text_align)
    latex.DrawLatex(pos_x, pos_y, CMS_TEXT)

    # extra text (preliminary, etc..)
    latex.SetTextFont(EXTRA_TEXT_FONT)
    latex.SetTextSize(EXTRA_TEXT_SIZE * PAD['t'])
    latex.SetTextAlign(text_align)
    latex.DrawLatex(pos_x, pos_y - REL_EXTRA_DY * CMS_TEXT_SIZE * PAD['t'], EXTRA_TEXT)

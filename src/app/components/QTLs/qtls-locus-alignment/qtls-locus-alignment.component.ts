import { Component, OnInit } from '@angular/core';
import { FormGroup, FormControl, Validators } from '@angular/forms';
import { QTLsResultsService } from '../../../services/qtls-results.service';
import { MatDialog } from '@angular/material';
import { QTLsLocusAlignmentBoxplotsComponent } from '../qtls-locus-alignment-boxplots/qtls-locus-alignment-boxplots.component';
import { environment } from '../../../../environments/environment' 

import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import * as PlotlyJS from 'plotly.js/dist/plotly.js';
import { PlotlyModule } from 'angular-plotly.js';

PlotlyModule.plotlyjs = PlotlyJS;


// import * as Plotly from '../../../../node_modules/plotly.js/dist/plotly.js';

declare let $: any;

export interface PopulationGroup {
  namecode: string;
  name: string;
  subPopulations: SubPopulation[];
}

export interface SubPopulation {
  value: string;
  viewValue: string;
}

export interface ReferenceGene {
  gene_id: string;
  gene_symbol: string;
}

export interface GeneVariants {
  gene_id: string;
  gene_symbol: string;
  rsnum: string;
}

@NgModule({
  imports: [CommonModule, PlotlyModule],
})

@Component({
  selector: 'app-qtls-locus-alignment',
  templateUrl: './qtls-locus-alignment.component.html',
  styleUrls: ['./qtls-locus-alignment.component.css']
})
export class QTLsLocusAlignmentComponent implements OnInit {

  locusAlignmentData: Object[];
  locusAlignmentDataR2: Object[];
  locusAlignmentDataR2NA: Object[]; 
  locusAlignmentDataRC: Object[];
  locusAlignmentDataQTopAnnot: Object;
  locusAlignmentGWASScatterData: Object[];
  locusAlignmentGWASScatterTitle: string;
  GWASData: Object[];
  GWASDataR2: Object[];
  GWASDataR2NA: Object[];
  selectedPop: string[];
  selectedGene: string;
  selectedGeneSymbol: string;
  selectedDist: string;
  public graph = null;
  public scatter = null;
  showPopover: boolean;
  collapseInput: boolean;
  popoverData: Object;
  disableLocusQuantification: boolean;
  requestID: number;
  associationFile: string;
  quantificationFile: string;
  genotypeFile: string;
  gwasFile: string;
  LDFile: string;
  recalculateAttempt: string;
  recalculatePopAttempt: string;
  recalculateGeneAttempt: string;
  recalculateDistAttempt: string;
  recalculateRefAttempt: string;
  newSelectedPop: string;
  newSelectedGene: string;
  newSelectedDist: string;
  blurLoadMain: boolean;
  disableInputs: boolean;
  selectedPvalThreshold: number;

  GWASScatterThreshold = new FormGroup({
    pvalThreshold: new FormControl("1.0", [Validators.pattern("^(\-?[0-9]*\.?[0-9]*)$"), Validators.min(0.0), Validators.max(1.0)])
  });

  select_qtls_samples: string;
  select_gwas_sample: string;

  constructor(private data: QTLsResultsService, public dialog: MatDialog) { }

  ngOnInit() {
    this.data.currentBlurLoadMain.subscribe(blurLoadMain => this.blurLoadMain = blurLoadMain);
    this.data.currentCollapseInput.subscribe(collapseInput => {
      this.collapseInput = collapseInput;
      // ensure graphs are properly positioned when data input panel collapse is toggled
      if (this.collapseInput) { // input panel collapsed
        if ($("#qtls-locus-alignment-plot").hasClass("justify-content-start")) {
          $("#qtls-locus-alignment-plot").addClass("justify-content-center");
          $("#qtls-locus-alignment-plot").removeClass("justify-content-start");
        }
        // if ($("#qtls-locus-alignment-scatter-plot").hasClass("justify-content-start")) {
        //   $("#qtls-locus-alignment-scatter-plot").addClass("justify-content-center");
        //   $("#qtls-locus-alignment-scatter-plot").removeClass("justify-content-start");
        // }
      } else { // input panel shown
        if ($("#qtls-locus-alignment-plot").hasClass("justify-content-center")) {
          $("#qtls-locus-alignment-plot").addClass("justify-content-start");
          $("#qtls-locus-alignment-plot").removeClass("justify-content-center");
        }
        // if ($("#qtls-locus-alignment-scatter-plot").hasClass("justify-content-center")) {
        //   $("#qtls-locus-alignment-scatter-plot").addClass("justify-content-start");
        //   $("#qtls-locus-alignment-scatter-plot").removeClass("justify-content-center");
        // }
      }
    });
    this.disableInputs = false;
    this.showPopover = false;

    this.data.currentLocusQuantification.subscribe(disableLocusQuantification => {
      this.disableLocusQuantification = disableLocusQuantification;
    });
    this.data.currentMainData.subscribe(mainData => {
      if (mainData) {
        this.recalculateAttempt = mainData["info"]["recalculateAttempt"][0]; // recalculation attempt ?
        this.recalculatePopAttempt = mainData["info"]["recalculatePop"][0]; // recalculation attempt when pop changed ?
        this.recalculateGeneAttempt = mainData["info"]["recalculateGene"][0]; // recalculation attempt when gene changed ?
        this.recalculateDistAttempt = mainData["info"]["recalculateDist"][0]; // recalculation attempt when cis-QTL distance changed ?
        this.recalculateRefAttempt = mainData["info"]["recalculateRef"][0]; // recalculation attempt when ref rsnum changed ?
        this.select_qtls_samples = mainData["info"]["select_qtls_samples"][0]; // use QTLs sample data files ?
        this.select_gwas_sample = mainData["info"]["select_gwas_sample"][0]; // use GWAS sample data file ?
        this.associationFile = mainData["info"]["inputs"]["association_file"][0]; // association filename
        this.quantificationFile = mainData["info"]["inputs"]["quantification_file"][0]; // quantification filename
        this.genotypeFile = mainData["info"]["inputs"]["genotype_file"][0]; // genotype filename
        this.gwasFile = mainData["info"]["inputs"]["gwas_file"][0] // gwas filename
        this.LDFile = mainData["info"]["inputs"]["ld_file"][0] // LD filename
        this.newSelectedPop = mainData["info"]["inputs"]["select_pop"][0]; // inputted populations
        this.newSelectedGene = mainData["info"]["inputs"]["select_gene"][0]; // inputted gene
        this.newSelectedDist = mainData["info"]["inputs"]["select_dist"][0]; // inputted cis-QTL distance
        this.requestID = mainData["info"]["inputs"]["request"][0]; // request id
        this.locusAlignmentData = mainData["locus_alignment"]["data"][0]; // locus alignment data
        this.locusAlignmentDataRC = mainData["locus_alignment"]["rc"][0]; // locus alignment RC data
        this.locusAlignmentDataQTopAnnot = mainData["locus_alignment"]["top"][0][0]; // locus alignment Top Gene data
        this.locusAlignmentGWASScatterData = mainData["locus_alignment_gwas_scatter"]["data"][0]; // locus alignment scatter data
        this.locusAlignmentGWASScatterTitle = mainData["locus_alignment_gwas_scatter"]["title"][0]; // locus alignment scatter title
        this.GWASData = mainData["gwas"]["data"][0]; // gwas data

        var newSelectedPopList = this.newSelectedPop.split('+');
        this.selectedPop = newSelectedPopList; // recalculated new population selection
        this.selectedGene = this.newSelectedGene; // recalculated new gene selection
        this.selectedDist = this.newSelectedDist;

        this.selectedGeneSymbol = this.locusAlignmentDataQTopAnnot["gene_symbol"];

        if (this.locusAlignmentData && this.locusAlignmentData[0]) {
          // differentiate variants with/without R2 data for popovers
          this.locusAlignmentDataR2 = this.getPopoverData(this.locusAlignmentData);
          this.locusAlignmentDataR2NA = this.getPopoverDataR2NA(this.locusAlignmentData);
          this.GWASDataR2 = this.getPopoverDataGWAS(this.GWASData);
          this.GWASDataR2NA = this.getPopoverDataGWASR2NA(this.GWASData);
          // check if there is data in GWAS object
          if ((this.GWASData && this.GWASData[0]) && (this.locusAlignmentGWASScatterData && this.locusAlignmentGWASScatterData[0])) {
            // if there is, graph GWAS plot and scatter plot
            this.locusAlignmentPlotGWAS(this.locusAlignmentDataR2, this.locusAlignmentDataR2NA, this.GWASDataR2, this.GWASDataR2NA, this.locusAlignmentDataRC, this.locusAlignmentDataQTopAnnot);
            this.selectedPvalThreshold = 1.0;
            this.locusAlignmentScatterPlot(this.locusAlignmentGWASScatterData, this.locusAlignmentGWASScatterTitle, this.selectedPvalThreshold);
          } else {
            // if not, do not graph GWAS plot or scatter plot
            this.locusAlignmentPlot(this.locusAlignmentDataR2, this.locusAlignmentDataR2NA, this.locusAlignmentDataRC, this.locusAlignmentDataQTopAnnot)
            // this.locusAlignmentPlot(this.locusAlignmentData, this.locusAlignmentDataRC, this.locusAlignmentDataQTopAnnot);
          }
        }

      }
    });
  }

  getPopoverData(geneData) {
    var dataR2 = [];
    for (var i = 0; i < geneData.length; i++) {
      if ("R2" in geneData[i] && geneData[i]["R2"] != "NA") {
        dataR2.push(geneData[i]);
      }
    }
    return dataR2;
  }

  getPopoverDataR2NA(geneData) {
    var dataR2NA = [];
    for (var i = 0; i < geneData.length; i++) {
      if (!("R2" in geneData[i]) || geneData[i]["R2"] == "NA") {
        dataR2NA.push(geneData[i]);
      }
    }
    return dataR2NA;
  }

  getPopoverDataGWAS(geneData) {
    var dataGWASR2 = [];
    for (var i = 0; i < geneData.length; i++) {
      if ("R2" in geneData[i] && geneData[i]["R2"] != "NA") {
        dataGWASR2.push(geneData[i]);
      }
    }
    return dataGWASR2;
  }

  getPopoverDataGWASR2NA(geneData) {
    var dataGWASR2NA = [];
    for (var i = 0; i < geneData.length; i++) {
      if (!("R2" in geneData[i]) || geneData[i]["R2"] == "NA") {
        dataGWASR2NA.push(geneData[i]);
      }
    }
    return dataGWASR2NA;
  }

  getXData(geneData) {
    var xData = [];
    for (var i = 0; i < geneData.length; i++) {
      xData.push(geneData[i]['pos'] / 1000000.0);
    }
    return xData;
  }

  getYData(geneData) {
    var yData = [];
    for (var i = 0; i < geneData.length; i++) {
      yData.push(Math.log10(geneData[i]['pval_nominal']) * -1.0);
    }
    return yData;
  }

  getHoverData(geneData) {
    var hoverData = [];
    for (var i = 0; i < geneData.length; i++) {
      if ('rsnum' in geneData[i]) {
        hoverData.push('chr' + geneData[i]['chr'] + ':' + geneData[i]['pos'] + '<br>' + geneData[i]['rsnum'] + '<br>' + 'Ref/Alt: ' + geneData[i]['ref'] + '/' + geneData[i]['alt'] + '<br>' + '<i>P</i>-value: ' + geneData[i]['pval_nominal'] + '<br>' + 'Slope: ' + geneData[i]['slope'] + '<br>' + "R2: " + (geneData[i]['R2'] ? geneData[i]['R2'] : "NA").toString());
      } else {
        hoverData.push('chr' + geneData[i]['chr'] + ':' + geneData[i]['pos'] + '<br>' + 'Ref/Alt: ' + geneData[i]['ref'] + '/' + geneData[i]['alt'] + '<br>' + '<i>P</i>-value: ' + geneData[i]['pval_nominal'] + '<br>' + 'Slope: ' + geneData[i]['slope'] + '<br>' + "R2: " + (geneData[i]['R2'] ? geneData[i]['R2'] : "NA").toString());
      }
    }
    return hoverData;
  }

  getHoverDataGWAS(geneGWASData) {
    var hoverData = [];
    for (var i = 0; i < geneGWASData.length; i++) {
      if ('rsnum' in geneGWASData[i]) {
        hoverData.push('chr' + geneGWASData[i]['chr'] + ':' + geneGWASData[i]['pos'] + '<br>' + geneGWASData[i]['rsnum'] + '<br>' + 'Ref/Alt: ' + geneGWASData[i]['ref'] + '/' + geneGWASData[i]['alt'] + '<br>' + '<i>P</i>-value: ' + geneGWASData[i]['pvalue'] + '<br>' + 'Slope: ' + geneGWASData[i]['slope'] + '<br>' + "R2: " + (geneGWASData[i]['R2'] ? geneGWASData[i]['R2'] : "NA").toString());
      } else {
        hoverData.push('chr' + geneGWASData[i]['chr'] + ':' + geneGWASData[i]['pos'] + '<br>' + 'Ref/Alt: ' + geneGWASData[i]['ref'] + '/' + geneGWASData[i]['alt'] + '<br>' + '<i>P</i>-value: ' + geneGWASData[i]['pvalue'] + '<br>' + 'Slope: ' + geneGWASData[i]['slope'] + '<br>' + "R2: " + (geneGWASData[i]['R2'] ? geneGWASData[i]['R2'] : "NA").toString());
      }
    }
    return hoverData;
  }

  getHoverDataRC(geneDataRC) {
    var hoverDataRC = [];
    for (var i = 0; i < geneDataRC.length; i++) {
      hoverDataRC.push('chr' + geneDataRC[i]['chr'] + ':' + geneDataRC[i]['pos'] + '<br>' + 'Rate: ' + geneDataRC[i]['rate']);
    }
    return hoverDataRC;
  }

  getYGWASData(geneData) {
    var yData = [];
    for (var i = 0; i < geneData.length; i++) {
      yData.push(Math.log10(geneData[i]['pvalue']) * -1.0);
    }
    return yData;
  }

  getYLDRefGWAS(xData, yGWASData) {
    var LDRefPos = this.locusAlignmentDataQTopAnnot['pos'] / 1000000.0;
    var pointIndex = xData.indexOf(LDRefPos);
    var yLDRef = yGWASData[pointIndex];
    return yLDRef;
  }

  getColorData(geneData) {
    var colorData = [];
    for (var i = 0; i < geneData.length; i++) {
      colorData.push(geneData[i]['R2']);
    }
    // normalize R2 color data between 0 and 1 for color spectrum
    for (var i = 0; i < colorData.length; i++) {
      colorData[i] = (colorData[i] / (Math.max.apply(null, colorData) / 100) ) / 100.0;
    }
    return colorData;
  }

  getXDataRC(geneDataRC) {
    var xData = [];
    for (var i = 0; i < geneDataRC.length; i++) {
      xData.push(geneDataRC[i]['pos'] / 1000000.0);
    }
    return xData;
  }

  getYDataRC(geneDataRC) {
    var yData = [];
    for (var i = 0; i < geneDataRC.length; i++) {
      yData.push(geneDataRC[i]['rate']);
    }
    return yData;
  }

  getMaxFinite(data) {
    var max = Number.NEGATIVE_INFINITY;
    for(var i = 0; 0 < data.length; i++) {
      if (isFinite(data[i])) {
        if (data[i] > max) {
          max = data[i];
        }
      }
    }
    return max;
  }

  locusAlignmentPlotGWAS(geneDataR2, geneDataR2NA, geneGWASDataR2, geneGWASDataR2NA, geneDataRC, qDataTopAnnot) {
    var xData = this.getXData(geneDataR2);
    var yData = this.getYData(geneDataR2);
    var xDataR2NA = this.getXData(geneDataR2NA);
    var yDataR2NA = this.getYData(geneDataR2NA);
    var yGWASData = this.getYGWASData(geneGWASDataR2);
    var yGWASDataR2NA = this.getYGWASData(geneGWASDataR2NA);
    var colorData = this.getColorData(geneDataR2);
    var xDataRC = this.getXDataRC(geneDataRC);
    var yDataRC = this.getYDataRC(geneDataRC);
    var hoverData = this.getHoverData(geneDataR2);
    var hoverDataR2NA = this.getHoverData(geneDataR2NA);
    var hoverDataGWAS = this.getHoverDataGWAS(geneGWASDataR2);
    var hoverDataGWASR2NA = this.getHoverDataGWAS(geneGWASDataR2NA);
    var hoverDataRC = this.getHoverDataRC(geneDataRC);
    var xLDRef = qDataTopAnnot['pos'] / 1000000.0;
    var yLDRef = Math.log10(qDataTopAnnot['pval_nominal']) * -1.0;
    var yDataFinites = this.removeInfinities(yData);
    var topPvalY = Math.max.apply(null, yDataFinites);
    var topPvalIdx = yData.indexOf(topPvalY);
    // var yGWASDataFinites = this.removeInfinities(yGWASData);
    // var topPvalYGWAS = Math.max.apply(null, yGWASDataFinites);
    // var topPvalIdxGWAS = yGWASData.indexOf(topPvalYGWAS);
    // mark point with most significant P-value
    var topPvalMarker = {
      x: [xData[topPvalIdx]], 
      y: [topPvalY], 
      hoverinfo: 'none', 
      mode: 'markers', 
      type: 'scatter', 
      marker: {
        symbol: "diamond",
        size: 15,
        opacity: 1.0,
        color: "#ff66cc",
        line: {
          color: '#ff66cc',
          width: 2
        },
      },
      yaxis: 'y3'
    };
    // mark GWAS point with most significant P-value
    // var topPvalMarkerGWAS = {
    //   x: [xData[topPvalIdxGWAS]],
    //   y: [topPvalYGWAS],
    //   hoverinfo: 'none',
    //   mode: 'markers',
    //   type: 'scatter',
    //   marker: {
    //     symbol: "diamond",
    //     size: 14,
    //     opacity: 1.0,
    //     color: "red",
    //     // color: [colorData[topIdx]],
    //     // colorscale: 'Viridis',
    //     // reversescale: true,
    //     line: {
    //       color: 'red',
    //       width: 2
    //     },
    //   }
    // };
    // highlight top point
    var LDRefHighlight = {
      x: [xLDRef],
      y: [yLDRef],
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 15,
        color: "red"
      },
      yaxis: 'y3'
    };
    // console.log("qDataTopAnnot", qDataTopAnnot);
    // console.log("xLDRef", xLDRef);
    // console.log("yLDRef", yLDRef);
    // console.log("LDRefHighlight", LDRefHighlight);
    // highlight top point GWAS
    var LDRefHighlightGWAS = {
      x: [xLDRef],
      y: [this.getYLDRefGWAS(xData, yGWASData)],
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 15,
        color: "red"
      },
      yaxis: 'y2'
    };
    // graph GWAS scatter where R2 != NA
    var trace1 = {
      x: xData,
      y: yGWASData,
      text: hoverDataGWAS,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        color: colorData,
        colorscale: 'Viridis',
        reversescale: true,
        line: {
          color: 'black',
          width: 1
        },
      },
      yaxis: 'y2'
    };
    // graph GWAS scatter where R2 == NA
    var trace1R2NA = {
      x: xDataR2NA,
      y: yGWASDataR2NA,
      text: hoverDataGWASR2NA,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        color: '#cccccc',
        line: {
          color: 'black',
          width: 1
        },
      },
      yaxis: 'y2'
    };
    // graph recombination rate line GWAS
    var trace2 = {
      x: xDataRC,
      y: yDataRC,
      text: hoverDataRC,
      hoverinfo: 'text',
      yaxis: 'y4',
      type: 'scatter',
      opacity: 0.35,
      line: {
        color: 'blue',
        width: 1
      }
    };
    // graph scatter where R2 != NA
    var trace3 = {
      x: xData,
      y: yData,
      text: hoverData,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        color: colorData,
        colorscale: 'Viridis',
        reversescale: true,
        line: {
          color: 'black',
          width: 1
        },
      },
      yaxis: 'y3'
    };
    // graph scatter where R2 == NA
    var trace3R2NA = {
      x: xDataR2NA,
      y: yDataR2NA,
      text: hoverDataR2NA,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        color: '#cccccc',
        line: {
          color: 'black',
          width: 1
        },
      },
      yaxis: 'y3'
    };
    // graph recombination rate line
    var trace4 = {
      x: xDataRC,
      y: yDataRC,
      text: hoverDataRC,
      hoverinfo: 'text',
      yaxis: 'y5',
      type: 'scatter',
      opacity: 0.35,
      line: {
        color: 'blue',
        width: 1
      }
    };
    // graph gene density where R2 != NA
    var trace5 = {
      x: xData,
      y: Array(xData.length).fill(0),
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        symbol: "line-ns-open",
        size: 16,
        color: colorData,
        colorscale: 'Viridis',
        reversescale: true
      },
      yaxis: 'y'
    };
    // graph gene density where R2 == NA
    var trace5R2NA = {
      x: xDataR2NA,
      y: Array(xDataR2NA.length).fill(0),
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        symbol: "line-ns-open",
        size: 16,
        color: '#cccccc',
      },
      yaxis: 'y'
    };
    var pdata = [topPvalMarker, LDRefHighlight, LDRefHighlightGWAS, trace1, trace1R2NA, trace2, trace3, trace3R2NA, trace4, trace5, trace5R2NA];
    // var pdata = [topPvalMarker, topPvalMarkerGWAS, LDRefHighlight, LDRefHighlightGWAS, trace1, trace2, trace3, trace4];
    var chromosome = qDataTopAnnot['chr'];
    var playout = {
      title: {
        text: 'QTLs-GWAS Chromosome ' + chromosome + ' Variants',
        xref: 'paper'
      },
      font: {
        color: 'black'
      },
      width: 1000,
      height: 1180,
      yaxis: {
        autorange: true,
        fixedrange: true,
        // overlaying: false,
        // title: "Gene Density",
        domain: [0, 0.025],
        zeroline: false,
        showgrid: false,
        showticklabels: false,
        linecolor: 'black',
        linewidth: 1,
        mirror: true,
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      yaxis2: {
        autorange: true,
        automargin: true,
        title: "GWAS -log10(<i>P</i>-value)",
        domain: [0.03, 0.54],
        zeroline: false,
        linecolor: 'black',
        linewidth: 1,
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      yaxis3: {
        autorange: true,
        automargin: true,
        title: "QTLs -log10(<i>P</i>-value), " + this.selectedGeneSymbol,
        domain: [0.56, 1],
        zeroline: false,
        linecolor: 'black',
        linewidth: 1,
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      yaxis4: {
        autorange: true,
        automargin: true,
        title: 'Recombination Rate (cM/Mb)',
        titlefont: {
          color: 'blue'
        },
        tickfont: {
          color: 'blue'
        },
        overlaying: 'y2',
        side: 'right',
        showgrid: false,
        zeroline: false,
        linecolor: 'black',
        linewidth: 1
      },
      yaxis5: {
        autorange: true,
        automargin: true,
        title: 'Recombination Rate (cM/Mb)',
        titlefont: {
          color: 'blue'
        },
        tickfont: {
          color: 'blue'
        },
        overlaying: 'y3',
        side: 'right',
        showgrid: false,
        zeroline: false,
        linecolor: 'black',
        linewidth: 1
      },
      xaxis: {
        autorange: true,
        title: "Chromosome " + chromosome + " (Mb)",
        zeroline: false,
        linecolor: 'black',
        linewidth: 1,
        mirror: "allticks",
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      images: [
        {
          x: 0,
          y: 1.01,
          sizex: 1.0,
          sizey: 1.0,
          source: environment.endpoint + "assets/images/qtls_locus_alignment_r2_legend_transparent.png",
          xanchor: "left",
          xref: "paper",
          yanchor: "bottom",
          yref: "paper"
        }
      ],
      margin: {
        l: 40,
        r: 40,
        b: 80,
        t: 120
      },
      showlegend: false,
      clickmode: 'event',
      hovermode: 'closest',
    };
    this.graph = {
      data: pdata,
      layout: playout, 
      config: {
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"],
        toImageButtonOptions: {
          format: 'svg', // one of png, svg, jpeg, webp
          filename: 'locus_alignment_manhattan_gwas',
          width: 1000,
          height: 1100,
          scale: 1 // Multiply title/legend/axis/canvas sizes by this factor
        }
      }
    };
  }

  locusAlignmentPlot(geneDataR2, geneDataR2NA, geneDataRC, qDataTopAnnot) {
    var xData = this.getXData(geneDataR2);
    var yData = this.getYData(geneDataR2);
    var xDataR2NA = this.getXData(geneDataR2NA);
    var yDataR2NA = this.getYData(geneDataR2NA);
    var colorData = this.getColorData(geneDataR2);
    var xDataRC = this.getXDataRC(geneDataRC);
    var yDataRC = this.getYDataRC(geneDataRC);
    var hoverData = this.getHoverData(geneDataR2);
    var hoverDataR2NA = this.getHoverData(geneDataR2NA);
    var hoverDataRC = this.getHoverDataRC(geneDataRC);
    var yDataFinites = this.removeInfinities(yData);
    var topPvalY = Math.max.apply(null, yDataFinites);
    var topPvalIdx = yData.indexOf(topPvalY);
    // mark point with most significant P-value
    var topPvalMarker = {
      x: [xData[topPvalIdx]],
      y: [topPvalY],
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        symbol: "diamond",
        size: 15,
        opacity: 1.0,
        color: "#ff66cc",
        line: {
          color: '#ff66cc',
          width: 2
        },
      },
      yaxis: 'y2'
    };
    // highlight top point
    var LDRefHighlight = {
      x: [qDataTopAnnot['pos'] / 1000000.0],
      y: [Math.log10(qDataTopAnnot['pval_nominal']) * -1.0],
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        opacity: 1.0, 
        size: 15,
        color: "red"
      },
      yaxis: 'y2'
    };

    // console.log("qDataTopAnnot", qDataTopAnnot);
    // console.log("xLDRef", qDataTopAnnot['pos'] / 1000000.0);
    // console.log("yLDRef", Math.log10(qDataTopAnnot['pval_nominal']) * -1.0);
    // console.log("LDRefHighlight", LDRefHighlight);

    // graph scatter where R2 != NA
    var trace1 = {
      x: xData,
      y: yData,
      text: hoverData,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        opacity: 1.0,
        color: colorData,
        colorscale: 'Viridis',
        reversescale: true,
        line: {
          color: 'black',
          width: 1
        },
      },
      yaxis: 'y2'
    };
    // graph scatter where R2 == NA
    var trace1R2NA = {
      x: xDataR2NA,
      y: yDataR2NA,
      text: hoverDataR2NA,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        opacity: 1.0,
        color: '#cccccc',
        line: {
          color: 'black',
          width: 1
        },
      },
      yaxis: 'y2'
    };
    // graph recombination rate line
    var trace2 = {
      x: xDataRC,
      y: yDataRC,
      text: hoverDataRC,
      hoverinfo: 'text',
      yaxis: 'y3',
      type: 'scatter',
      opacity: 0.35,
      line: {
        color: 'blue',
        width: 1
      }
    };
    // graph gene density where R2 != NA
    var trace3 = {
      x: xData,
      y: Array(xData.length).fill(0),
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        symbol: "line-ns-open",
        size: 16,
        color: colorData,
        colorscale: 'Viridis',
        reversescale: true
      },
      yaxis: 'y'

    };
    // graph gene density where R2 == NA
    var trace3R2NA = {
      x: xDataR2NA,
      y: Array(xDataR2NA.length).fill(0),
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        symbol: "line-ns-open",
        size: 16,
        color: '#cccccc',
      },
      yaxis: 'y'

    };
    var pdata = [topPvalMarker, LDRefHighlight, trace1, trace1R2NA, trace2, trace3, trace3R2NA];
    // round most significant pval to next whole number
    var chromosome = qDataTopAnnot['chr'];
    var playout = {
      title: {
        text: 'QTLs Chromosome ' + chromosome + ' Variants',
        xref: 'paper'
      },
      width: 1000,
      height: 780,
      yaxis: {
        autorange: true,
        fixedrange: true,
        // overlaying: false,
        // title: "Gene Density",
        domain: [0, 0.04],
        zeroline: false,
        showgrid: false,
        showticklabels: false,
        linecolor: 'black',
        linewidth: 1,
        mirror: true,
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      yaxis2: {
        autorange: true,
        automargin: true,
        // overlaying: 'y3',
        title: "QTLs -log10(<i>P</i>-value), " + this.selectedGeneSymbol,
        domain: [0.05, 1],
        zeroline: false,
        linecolor: 'black',
        linewidth: 1,
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      yaxis3: {
        autorange: true,
        automargin: true,
        overlaying: 'y2',
        title: 'QTLs Recombination Rate (cM/Mb)',
        titlefont: {
          color: 'blue'
        },
        tickfont: {
          color: 'blue'
        },
        side: 'right',
        showgrid: false,
        zeroline: false,
        linecolor: 'black',
        linewidth: 1
      },
      xaxis: {
        autorange: true,
        title: "Chromosome " + chromosome + " (Mb)",
        // dtick: 0.1,
        zeroline: false,
        linecolor: 'black',
        linewidth: 1,
        mirror: "allticks",
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      images: [
        {
          x: 0,
          y: 1.02,
          sizex: 1.0,
          sizey: 1.0,
          source: environment.endpoint + "assets/images/qtls_locus_alignment_r2_legend_transparent.png",
          xanchor: "left",
          xref: "paper",
          yanchor: "bottom",
          yref: "paper"
        }
      ],
      margin: {
        l: 40,
        r: 40,
        b: 80,
        t: 120
      },
      showlegend: false,
      clickmode: 'event',
      hovermode: 'closest'
    };
    this.graph = {
      data: pdata,
      layout: playout,
      config: {
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"],
        toImageButtonOptions: {
          format: 'svg', // one of png, svg, jpeg, webp
          filename: 'locus_alignment_manhattan',
          width: 1000,
          height: 700,
          scale: 1 // Multiply title/legend/axis/canvas sizes by this factor
        }
      }
    };
    // var pconfig = {
    //   displaylogo: false,
    //   modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"],
    //   toImageButtonOptions: {
    //     format: 'svg', // one of png, svg, jpeg, webp
    //     filename: 'locus_alignment_manhattan',
    //     width: 1000,
    //     height: 700,
    //     scale: 1 // Multiply title/legend/axis/canvas sizes by this factor
    //   }
    // };
    // PlotlyJS.newPlot("qtls-locus-alignment-plot", pdata, playout, pconfig);
    // PlotlyJS.moveTraces("qtls-locus-alignment-plot", 1, 0);
    // PlotlyJS.moveTraces("qtls-locus-alignment-plot", 2, 0);
    // var manhattanPlot = $('#qtls-locus-alignment-plot')[0];
    // // var helper = this.clickPoint(event);
    // manhattanPlot.on('plotly_click', function(data){
    //   this.clickPoint(data);
    //   var dat = 0;
    // });
  }

  closePopover() {
    $('.popover').hide();
    this.showPopover = false;
  }
  
  // close popover if click anywhere else
  closePopover2(event) {
    if (event.points == null && this.showPopover == true) {
      this.showPopover = false;
    } else {
      this.showPopover = false;
      $('.popover').hide();
    }
  }

  async makeLDRef() {
    var selectedGeneString = this.selectedGene;
    var selectedPopString = this.selectedPop.join('+');
    var selectedDistNumber = this.selectedDist;
    var selectedRefString = this.popoverData["rsnum"];
    var recalculateAttempt = "true";
    var recalculatePop = "false";
    var recalculateGene = "false";
    var recalculateDist = "false";
    var recalculateRef = "true";
    // reset
    this.closePopover();
    this.data.changeBlurLoadMain(true);
    this.data.changeECAVIARData(null);
    this.data.changeHyprcolocData(null);
    this.disableInputs = true;
    $("#ldref-search-warning").hide();
    $(".blur-loading-main").addClass("blur-overlay");
    $(".blur-loading-ecaviar").addClass("blur-overlay");
    // calculate
    this.data.recalculateMain(this.select_qtls_samples, this.select_gwas_sample, this.associationFile, this.quantificationFile, this.genotypeFile, this.gwasFile, this.LDFile, this.requestID, selectedPopString, selectedGeneString, selectedDistNumber, selectedRefString, recalculateAttempt, recalculatePop, recalculateGene, recalculateDist, recalculateRef)
      .subscribe(
        res => {
          this.data.changeMainData(res);
          this.data.changeBlurLoadMain(false);
          this.disableInputs = false;
          $(".blur-loading-main").removeClass("blur-overlay");
          this.recalculatePopAttempt = "false";
          this.recalculateGeneAttempt = "false";
          this.recalculateDistAttempt = "false";
          this.recalculateRefAttempt = "false";
        },
        error => {
          this.handleError(error);
          this.data.changeBlurLoadMain(false);
          this.disableInputs = false;
          $(".blur-loading-main").removeClass("blur-overlay");
        }
      );
  }

  linkLDpop() {
    var selectedRefString = this.popoverData["rsnum"];
    var QTopAnnotRef = this.locusAlignmentDataQTopAnnot["rsnum"];
    var selectedPopString = this.selectedPop.join('%2B');
    var url = "https://ldlink.nci.nih.gov/?tab=ldpop&var1=" + selectedRefString + "&var2=" + QTopAnnotRef + "&pop=" + selectedPopString + "&r2_d=r2"
    var win = window.open(url, '_blank');
    win.focus();
  } 

  linkGWAS() {
    var selectedRefString = this.popoverData["rsnum"];
    var url = "https://www.ebi.ac.uk/gwas/search?query=" + selectedRefString
    var win = window.open(url, '_blank');
    win.focus();
  }

  linkGnomADBrowser() {
    var chromosome = this.popoverData["chr"];
    var position = this.popoverData["pos"];
    var ref = this.popoverData["ref"];
    var alt = this.popoverData["alt"];
    var url = "http://gnomad.broadinstitute.org/variant/" + chromosome + "-" + position + "-" + ref + "-" + alt
    var win = window.open(url, '_blank');
    win.focus();
  }

  async showBoxplot() {
    if (this.popoverData) {
      // console.log("LOCUS ALIGNMENT BOXPLOT DIALOG MODULE OPENED");
      this.dialog.open(QTLsLocusAlignmentBoxplotsComponent, {
        data: this.popoverData
      });
      this.closePopover();
      this.popoverData = null;
    }
  }

  clickPoint(event) {    
    if (event.points) {
      // console.log("curveNumber", event.points[0].curveNumber);
      // console.log("pointIndex", event.points[0].pointIndex);
      // console.log(event.points[0]);
      // if (event.points[0].hasOwnProperty("marker.color")) {
        // only show popovers for scatter points not recomb line (points w/ markers)
        var top = event.event.pointerY;
        var left = event.event.pointerX;
        // this.popoverData = this.locusAlignmentData[event.points[0].pointIndex];
        // console.log("popoverData", this.popoverData);
        if (this.GWASData[0]) { // if GWAS data is displayed in Manhattan plot
          // console.log("GWAS data displayed");
          if (event.points[0].curveNumber == 3 || event.points[0].curveNumber == 4) { // if GWAS data is clicked
            // console.log("GWAS data clicked.");
            if (event.points[0].curveNumber == 3) {
              this.popoverData = this.GWASDataR2[event.points[0].pointIndex];
            } else {
              this.popoverData = this.GWASDataR2NA[event.points[0].pointIndex];
            }
            $('.popover').show();
            if (this.collapseInput) { // input panel collapsed
              $('.popover').show().css({
                position: "absolute",
                top: top + 70, 
                left: left + 245
              });
            } else { // input panel shown
              $('.popover').show().css({
                position: "absolute",
                top: top + 70, 
                left: left + 25
              });
            }
            this.showPopover = true;
            // console.log("popoverData", this.popoverData);
          } else { // if Association data is clicked
            if (event.points[0].curveNumber == 6 || event.points[0].curveNumber == 7) {
              // console.log("Association data clicked.");
              if (event.points[0].curveNumber == 6) {
                var associationData = this.locusAlignmentDataR2[event.points[0].pointIndex];
              } else {
                var associationData = this.locusAlignmentDataR2NA[event.points[0].pointIndex];
              }
              this.popoverData = {
                chr: associationData["chr"], 
                pos: associationData["pos"], 
                variant_id: associationData["variant_id"], 
                gene_id: associationData["gene_id"], 
                gene_symbol: associationData["gene_symbol"],
                ref: associationData["ref"], 
                alt: associationData["alt"], 
                rsnum: associationData["rsnum"], 
                pvalue: associationData["pval_nominal"], 
                zscore: associationData["zscore"], 
                effect: associationData["effect"],  
                slope: associationData["slope"], 
                se: associationData["se"], 
                R2: associationData["R2"], 
                tss_distance: associationData["tss_distance"]
              };
              $('.popover').show();
              if (this.collapseInput) { // input panel collapsed
                $('.popover').show().css({
                  position: "absolute",
                  top: top + 70, 
                  left: left + 245
                });
              } else { // input panel shown
                $('.popover').show().css({
                  position: "absolute",
                  top: top + 70, 
                  left: left + 25
                });
              }
              this.showPopover = true;
              // console.log("popoverData", this.popoverData);
            }
          }
        } else { // if no GWAS data is disaplyed in Manhattan plot
          // console.log("No GWAS data displayed");
          if (event.points[0].curveNumber == 2 || event.points[0].curveNumber == 3) {
            if (event.points[0].curveNumber == 2) {
              var associationData = this.locusAlignmentDataR2[event.points[0].pointIndex];
            } else {
              var associationData = this.locusAlignmentDataR2NA[event.points[0].pointIndex];
            }
            this.popoverData = {
              chr: associationData["chr"], 
              pos: associationData["pos"], 
              variant_id: associationData["variant_id"], 
              gene_id: associationData["gene_id"], 
              gene_symbol: associationData["gene_symbol"],
              ref: associationData["ref"], 
              alt: associationData["alt"], 
              rsnum: associationData["rsnum"], 
              pvalue: associationData["pval_nominal"], 
              zscore: associationData["zscore"], 
              effect: associationData["effect"],  
              slope: associationData["slope"], 
              se: associationData["se"], 
              R2: associationData["R2"], 
              tss_distance: associationData["tss_distance"]
            };
            $('.popover').show();
            if (this.collapseInput) { // input panel collapsed
              $('.popover').show().css({
                position: "absolute",
                top: top + 70, 
                left: left + 245
              });
            } else { // input panel shown
              $('.popover').show().css({
                position: "absolute",
                top: top + 70, 
                left: left + 25
              });
            }
            this.showPopover = true;
          }
          // console.log("popoverData", this.popoverData);
        }
        // $('.popover').show();
        // if (this.collapseInput) { // input panel collapsed
        //   $('.popover').show().css({
        //     position: "absolute",
        //     top: top - 125, 
        //     left: left + 245
        //   });
        // } else { // input panel shown
        //   $('.popover').show().css({
        //     position: "absolute",
        //     top: top - 125, 
        //     left: left + 25
        //   });
        // }
        // this.showPopover = true;
      // }
    } else {
      this.closePopover2(event);
    }
  }
  
  handleError(error) {
    console.log(error);
    this.closePopover();
    var errorTrimmed = error.error.trim().split('\n');
    // var errorMessage = errorTrimmed.slice(1, errorTrimmed.length - 1).join(' ');
    var errorMessage = errorTrimmed[2];
    console.log(errorMessage);
    this.data.changeErrorMessage(errorMessage);
  }

  getScatterX(scatterData, threshold) {
    var p_values = [];
    for (var i = 0; i < scatterData.length; i++) {
      if ("R2" in scatterData[i] && scatterData[i]["R2"] != "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          p_values.push(Math.log10(scatterData[i]['pvalue']) * -1.0);
        }
      }
    }
    return p_values;
  }

  getScatterXR2NA(scatterData, threshold) {
    var p_values = [];
    for (var i = 0; i < scatterData.length; i++) {
      if (!("R2" in scatterData[i]) || scatterData[i]["R2"] == "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          p_values.push(Math.log10(scatterData[i]['pvalue']) * -1.0);
        }
      }
    }
    return p_values;
  }

  getScatterY(scatterData, threshold) {
    var pval_nominals = [];
    for (var i = 0; i < scatterData.length; i++) {
      if ("R2" in scatterData[i] && scatterData[i]["R2"] != "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          pval_nominals.push(Math.log10(scatterData[i]['pval_nominal']) * -1.0);
        }
      }
    }
    return pval_nominals;
  }

  getScatterYR2NA(scatterData, threshold) {
    var pval_nominals = [];
    for (var i = 0; i < scatterData.length; i++) {
      if (!("R2" in scatterData[i]) || scatterData[i]["R2"] == "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          pval_nominals.push(Math.log10(scatterData[i]['pval_nominal']) * -1.0);
        }
      }
    }
    return pval_nominals;
  }

  getScatterHoverData(scatterData, threshold) {
    var hoverData = [];
    for (var i = 0; i < scatterData.length; i++) {
      if ("R2" in scatterData[i] && scatterData[i]["R2"] != "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          if ('rs' in scatterData[i]) {
            hoverData.push('chr' + scatterData[i]['chr'] + ':' + scatterData[i]['pos'] + '<br>' + scatterData[i]['rs'] + '<br>' + '<i>P</i>-value: ' + scatterData[i]['pval_nominal'] + '<br>' + "R2: " + (scatterData[i]['R2'] ? scatterData[i]['R2'] : "NA").toString());
          } else {
            hoverData.push('chr' + scatterData[i]['chr'] + ':' + scatterData[i]['pos'] + '<br>' + '<i>P</i>-value: ' + scatterData[i]['pval_nominal'] + '<br>' + "R2: " + (scatterData[i]['R2'] ? scatterData[i]['R2'] : "NA").toString());
          }
        }
      }
    }
    return hoverData;
  }

  getScatterHoverDataR2NA(scatterData, threshold) {
    var hoverData = [];
    for (var i = 0; i < scatterData.length; i++) {
      if (!("R2" in scatterData[i]) || scatterData[i]["R2"] == "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          if ('rs' in scatterData[i]) {
            hoverData.push('chr' + scatterData[i]['chr'] + ':' + scatterData[i]['pos'] + '<br>' + scatterData[i]['rs'] + '<br>' + '<i>P</i>-value: ' + scatterData[i]['pval_nominal'] + '<br>' + "R2: " + (scatterData[i]['R2'] ? scatterData[i]['R2'] : "NA").toString());
          } else {
            hoverData.push('chr' + scatterData[i]['chr'] + ':' + scatterData[i]['pos'] + '<br>' + '<i>P</i>-value: ' + scatterData[i]['pval_nominal'] + '<br>' + "R2: " + (scatterData[i]['R2'] ? scatterData[i]['R2'] : "NA").toString());
          }
        }
      }
    }
    return hoverData;
  }

  getScatterColorData(scatterData, threshold) {
    var colorData = [];
    for (var i = 0; i < scatterData.length; i++) {
      if ("R2" in scatterData[i] && scatterData[i]["R2"] != "NA") {
        if (scatterData[i]['pvalue'] <= parseFloat(threshold) && scatterData[i]['pval_nominal'] <= parseFloat(threshold)) {
          colorData.push(scatterData[i]['R2']);
        }
      }
    }
    // normalize R2 color data between 0 and 1 for color spectrum
    for (var i = 0; i < colorData.length; i++) {
      colorData[i] = (colorData[i] / (Math.max.apply(null, colorData) / 100) ) / 100.0;
    }
    return colorData;
  }

  getSum(arr) {
    var sum = 0;
    for (var i = 0; i < arr.length; i++) {
      if (isFinite(arr[i])) {
        sum += arr[i];
      } else {
        // console.log("NOT FINITE", arr[i]);
      }
    }
    return sum;
  }

  removeInfinities(arr) {
    var finites = [];
    for (var i = 0; i < arr.length; i ++) {
      if (isFinite(arr[i])) {
        finites.push(arr[i]);
      } else {

      }
    }
    return finites;
  }

  getLeastSquaresLine(valuesXRaw, valuesYRaw) {
    if (valuesXRaw.length === 0 || valuesYRaw.length === 0 || (valuesXRaw.length != valuesYRaw.length)) {
      return [ [], [] ];
    }
    var valuesXY = valuesXRaw.map(function(e, i) {
      return [e, valuesYRaw[i]];
    });
    var valuesXYRemoveInfinities = [];
    valuesXY.map(function(e) {
      if (isFinite(e[0]) && isFinite(e[1])) {
        valuesXYRemoveInfinities.push(e);
      }
    });
    var valuesX = valuesXYRemoveInfinities.map(function(e) {
      return e[0];
    });
    var valuesY = valuesXYRemoveInfinities.map(function(e) {
      return e[1];
    });
    var sumX = 0;
    var sumY = 0;
    var sumXY = 0;
    var sumXX = 0;
    var count = 0;
    var x = 0;
    var y = 0;
    for (var i = 0; i < valuesX.length; i++) {
        x = valuesX[i];
        y = valuesY[i];
        sumX += x;
        sumY += y;
        sumXX += x * x;
        sumXY += x * y;
        count++;
    }
    var m = (count * sumXY - sumX * sumY) / (count * sumXX - sumX * sumX);
    var b = (sumY / count) - (m * sumX) / count;
    var resultX = [];
    var resultY = [];
    for (var j = 0; j < valuesX.length; j++) {
        x = valuesX[j];
        y = (x * m) + b;
        resultX.push(x);
        resultY.push(y);
    }
    return [resultX, resultY];
}

  getSumDiffAbsSquared(arr1, arr2) {
    var d = [];
    for (var i = 0; i < arr1.length; i++) {
      d.push(Math.pow(Math.abs(arr1[i] - arr2[i]), 2));
    }
    var sum = this.getSum(d);
    return sum;
  }

  multiplyArrays(x, y) {
    var xy = [];
    for (var i = 0; i < x.length; i++) {
      xy.push(x[i] * y[i]);
    }
    return xy;
  }

  recalculateSpearmanCorrelationTitle(xData, yData) {
    if (xData.length > 0 && yData.length > 0) {
      var sortedX = xData.slice().sort(function(a, b){ return b - a })
      var xRank = xData.slice().map(function(p){ return sortedX.indexOf(p) + 1 });
      var sortedY = yData.slice().sort(function(a, b){ return b - a })
      var yRank = yData.slice().map(function(p){ return sortedY.indexOf(p) + 1 });
      var sumSquaredDiffRanks = this.getSumDiffAbsSquared(xRank, yRank);
      var numer = 6.0 * sumSquaredDiffRanks;
      var denom = xData.length * (Math.pow(xData.length, 2) - 1)
      var rho = 1 - (numer / denom);
      return "Spearman rho=" + rho.toFixed(3);
    } else {
      return "Spearman rho=NA"
    }
  }

  recalculatePearsonCorrelationTitle(xData, yData) {
    if (xData.length > 0 && yData.length > 0) {
      var xDataFinites = this.removeInfinities(xData);
      var yDataFinites = this.removeInfinities(yData);
      var xMean = xDataFinites.reduce(function(a, b){ return a + b }) / xDataFinites.length;
      var yMean = yDataFinites.reduce(function(a, b){ return a + b }) / yDataFinites.length;
      var xDataMinusMean = xData.map(function(i){ return i - xMean });
      var yDataMinusMean = yData.map(function(i){ return i - yMean });
      var xy = this.multiplyArrays(xDataMinusMean, yDataMinusMean);
      var xDataMinusMeanSquared = xDataMinusMean.map(function(i){ return Math.pow(i, 2) })
      var yDataMinusMeanSquared = yDataMinusMean.map(function(i){ return Math.pow(i, 2) })
      var numer = this.getSum(xy);
      var xSumDataMinusMeanSquared = this.getSum(xDataMinusMeanSquared);
      var ySumDataMinusMeanSquared = this.getSum(yDataMinusMeanSquared);
      var denom = Math.sqrt(xSumDataMinusMeanSquared * ySumDataMinusMeanSquared);
      var r = numer / denom;
      return "Pearson's r=" + r.toFixed(3);
    } else {
      return "Pearson's r=NA";
    }
  }

  locusAlignmentScatterPlot(scatterData, scatterTitle, threshold) {
    // console.log("scatterData", scatterData);
    var xData = this.getScatterX(scatterData, threshold);
    // console.log("xData", xData);
    var yData = this.getScatterY(scatterData, threshold);
    // console.log("yData", yData);
    var xDataR2NA = this.getScatterXR2NA(scatterData, threshold);
    // console.log("xDataR2NA", xDataR2NA);
    var yDataR2NA = this.getScatterYR2NA(scatterData, threshold);
    // console.log("yDataR2NA", yDataR2NA);
    var scatterColorData = this.getScatterColorData(scatterData, threshold);
    // console.log("scatterColorData", scatterColorData);
    var hoverData = this.getScatterHoverData(scatterData, threshold);
    // console.log("hoverData", hoverData);
    var hoverDataR2NA = this.getScatterHoverDataR2NA(scatterData, threshold);
    // console.log("hoverDataR2NA", hoverDataR2NA);
    var trace1 = {
      x: xData,
      y: yData,
      text: hoverData,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        color: scatterColorData,
        colorscale: 'Viridis',
        reversescale: true,
        line: {
          color: 'black',
          width: 1
        },
      }
    };
    var trace1R2NA = {
      x: xDataR2NA,
      y: yDataR2NA,
      text: hoverDataR2NA,
      hoverinfo: 'text',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 7,
        color: '#cccccc',
        line: {
          color: 'black',
          width: 1
        },
      }
    };
    var least_sqaures = this.getLeastSquaresLine(xData, yData);

    var trace2 = {
      x: least_sqaures[0],
      y: least_sqaures[1],
      // hoverinfo: 'x+y',
      mode: 'lines',
      type: 'scatter',
      name: "lines",
      line: {
        color: "blue",
      }
    };

    var pdata = [
      trace1, 
      trace1R2NA, 
      trace2
    ];
    // var pdata = [trace1];
    var playout = {
      title: {
        text: "QTL-GWAS <i>P</i>-value Correlation: " + ((scatterTitle == "RECALCULATE") ? this.recalculateSpearmanCorrelationTitle(xData, yData) + ", " + this.recalculatePearsonCorrelationTitle(xData, yData) : "Spearman " + scatterTitle.split(', ')[0] + ", Pearson's " + scatterTitle.split(', ')[1]),
        xref: 'paper'
      },
      font: {
        color: 'black'
      },
      width: 700,
      height: 700,
      yaxis: {
        autorange: true,
        automargin: true,
        title: "-log10(QTL <i>P</i>-value), " + this.selectedGeneSymbol,
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      xaxis: {
        autorange: true,
        automargin: true,
        title: "-log10(GWAS <i>P</i>-value)",
        font: {
          color: 'black'
        },
        tickfont: {
          color: 'black'
        }
      },
      // margin: {
      //   l: 40,
      //   r: 40,
      //   b: 40,
      //   t: 40
      // },
      showlegend: false,
      clickmode: 'none',
      hovermode: 'closest'
    };
    this.scatter = {
      data: pdata,
      layout: playout,
      config: {
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"],
        toImageButtonOptions: {
          format: 'svg', // one of png, svg, jpeg, webp
          filename: 'locus_alignment_gwas_scatter',
          width: 1000,
          height: 700,
          scale: 1 // Multiply title/legend/axis/canvas sizes by this factor
        }
      }
    };
  }

  changePvalThreshold(event: any) {
    var threshold = event.target.value;
    if (threshold == 1.0) { 
      this.locusAlignmentScatterPlot(this.locusAlignmentGWASScatterData, this.locusAlignmentGWASScatterTitle, this.selectedPvalThreshold);
    } else {
      if (threshold >= 0.0 && threshold <= 1.0) {
        if (threshold.length > 0) {
          this.selectedPvalThreshold = threshold;
          this.locusAlignmentScatterPlot(this.locusAlignmentGWASScatterData, "RECALCULATE", this.selectedPvalThreshold);
        } else {
          this.locusAlignmentScatterPlot(this.locusAlignmentGWASScatterData, "RECALCULATE", 0.0);
        }
      } 
    }
  }

  clearPvalThreshold() {
    this.selectedPvalThreshold = 1.0;
    this.GWASScatterThreshold.value.pvalThreshold = '1.0';
    this.locusAlignmentScatterPlot(this.locusAlignmentGWASScatterData, this.locusAlignmentGWASScatterTitle, this.selectedPvalThreshold);
  }

  pvalThresholdErrorMsg() {
    var msg = "";
    if (this.GWASScatterThreshold.value.pvalThreshold > 1.0) {
      msg = "Threshold must be <= 1.0";
    } else if (this.GWASScatterThreshold.value.pvalThreshold < 0.0) {
      msg = "Threshold must be >= 0.0";
    } else {
      msg = "Invalid threshold";
    }
    return msg;
  }

}

import { Component, OnInit, Inject } from '@angular/core';
import { EqtlResultsService } from '../../services/eqtl-results.service';
import { PlotComponent } from 'angular-plotly.js';
import { MatDialog } from '@angular/material';
import { EqtlResultsLocuszoomBoxplotsComponent } from '../eqtl-results-locuszoom-boxplots/eqtl-results-locuszoom-boxplots.component';
import { environment } from '../../../environments/environment' 
import { FormControl, Validators } from '@angular/forms';

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


@Component({
  selector: 'app-eqtl-results-locuszoom',
  templateUrl: './eqtl-results-locuszoom.component.html',
  styleUrls: ['./eqtl-results-locuszoom.component.css']
})
export class EqtlResultsLocuszoomComponent implements OnInit {

  locuszoomData: Object;
  locuszoomDataRC: Object;
  locuszoomDataQTopAnnot: Object;
  locuszoomScatterData: Object;
  locuszoomScatterTitle: string;
  GWASData: Object;
  geneList: ReferenceGene[];
  topGeneVariants: GeneVariants[];
  allGeneVariants: GeneVariants[];
  allGeneVariantsOrganized: Object;
  selectedGene: string;
  populationGroups: PopulationGroup[];
  selectedPop: string[];
  selectedPopFinal: string[];
  public graph = null;
  public scatter = null;
  showPopover: boolean;
  collapseInput: boolean;
  selectedRef: string;
  populationSelectedAll: boolean;
  popoverData: Object;
  disableGeneExpressions: boolean;
  inputChanged: boolean;
  requestID: number;
  associationFile: string;
  expressionFile: string;
  genotypeFile: string;
  gwasFile: string;
  recalculateAttempt: string;
  recalculatePopAttempt: string;
  recalculateGeneAttempt: string;
  recalculateRefAttempt: string;
  newSelectedPop: string;
  newSelectedGene: string;
  newSelectedRef: string;
  rsnumSearch: string;
  warningMessage: string;

  rsnumber = new FormControl('', [Validators.pattern("^(rs[0-9]+)?$")]);


  constructor(private data: EqtlResultsService, public dialog: MatDialog) { }

  ngOnInit() {
    this.data.currentCollapseInput.subscribe(collapseInput => this.collapseInput = collapseInput);
    this.populationGroups = this.populatePopulationDropdown();
    this.selectedPopFinal = [];
    this.populationSelectedAll = false;
    this.inputChanged = false;
    this.rsnumSearch = "";
    this.warningMessage = "";

    this.data.currentGeneExpressions.subscribe(disableGeneExpressions => {
      this.disableGeneExpressions = disableGeneExpressions;
    });
    this.data.currentMainData.subscribe(mainData => {
      if (mainData) {
        this.recalculateAttempt = mainData["info"]["recalculateAttempt"][0]; // recalculation attempt ?
        this.recalculatePopAttempt = mainData["info"]["recalculatePop"][0]; // recalculation attempt when pop changed ?
        this.recalculateGeneAttempt = mainData["info"]["recalculateGene"][0]; // recalculation attempt when gene changed ?
        this.recalculateRefAttempt = mainData["info"]["recalculateRef"][0]; // recalculation attempt when ref rsnum changed ?
        this.associationFile = mainData["info"]["inputs"]["association_file"][0]; // association filename
        this.expressionFile = mainData["info"]["inputs"]["expression_file"][0]; // expression filename
        this.genotypeFile = mainData["info"]["inputs"]["genotype_file"][0]; // genotype filename
        this.gwasFile = mainData["info"]["inputs"]["gwas_file"][0] // gwas filename
        this.newSelectedPop = mainData["info"]["inputs"]["select_pop"][0]; // inputted populations
        this.newSelectedGene = mainData["info"]["inputs"]["select_gene"][0]; // inputted gene
        this.newSelectedRef = mainData["info"]["inputs"]["select_ref"][0]; // inputted ref
        this.requestID = mainData["info"]["inputs"]["request"][0]; // request id
        this.geneList = mainData["info"]["gene_list"]["data"][0]; // get gene list & populate ref gene dropdown
        this.allGeneVariants = mainData["info"]["all_gene_variants"]["data"][0]; // get list of all rsnums for all genes
        this.topGeneVariants = mainData["info"]["top_gene_variants"]["data"][0]; // get list of top rsnum for all genes
        this.locuszoomData = mainData["locuszoom"]["data"][0]; // locuszoom data
        this.locuszoomDataRC = mainData["locuszoom"]["rc"][0]; // locuszoom RC data
        this.locuszoomDataQTopAnnot = mainData["locuszoom"]["top"][0][0]; // locuszoom Top Gene data
        this.locuszoomScatterData = mainData["locuszoom_scatter"]["data"][0]; // locuszoom scatter data
        this.locuszoomScatterTitle = mainData["locuszoom_scatter"]["title"][0]; // locuszoom scatter title
        this.GWASData = mainData["gwas"]["data"][0]; // gwas data
      }
      if (this.locuszoomData) {
        // this.geneVariantList = this.populateGeneVariantList(this.locuszoomData);
        // check if there is data in GWAS object
        if (this.GWASData[0] && this.locuszoomScatterData[0]) {
          // if there is, graph GWAS plot and scatter plot
          this.graph = this.locuszoomPlotGWAS(this.locuszoomData, this.GWASData, this.locuszoomDataRC, this.locuszoomDataQTopAnnot);
          this.scatter = this.locuszoomScatterPlot(this.locuszoomScatterData, this.locuszoomScatterTitle);
        } else {
          // if not, do not graph GWAS plot or scatter plot
          this.graph = this.locuszoomPlot(this.locuszoomData, this.locuszoomDataRC, this.locuszoomDataQTopAnnot)
        }
      }
    });
    if (this.recalculatePopAttempt == "false") {
      this.selectedPop = ["CEU", "TSI", "FIN", "GBR", "IBS"]; // default population EUR
      this.returnPopulationGroupFinal();
    } else {
      var newSelectedPopList = this.newSelectedPop.split('+');
      this.selectedPop = newSelectedPopList; // recalculated new population selection
      // this.selectedGene = this.newSelectedGene; // recalculated new gene selection
      this.recalculatePopAttempt = "false";
      this.returnPopulationGroupFinal();
    }
    if (this.geneList) {
      if(this.recalculateGeneAttempt == "false") {
        this.selectedGene = this.locuszoomDataQTopAnnot["gene_id"]; // default reference gene
      } else {
        this.selectedGene = this.newSelectedGene; // recalculated new gene selection
        this.recalculateGeneAttempt = "false";
      }
    }
    if (this.recalculateRefAttempt == "false") {
      this.selectedRef = "false"; // default ref rsnum
      this.rsnumSearch = this.locuszoomDataQTopAnnot["rsnum"];
    } else {
      this.selectedRef = this.newSelectedRef; // recalculated new gene selection
      this.recalculateRefAttempt = "false";
      this.rsnumSearch = this.selectedRef;
    }
    if (this.allGeneVariants && this.geneList) {
      this.populateAllGeneVariantLists(this.allGeneVariants, this.geneList); // organize all eQTL variants by gene
    }
  }

  populateAllGeneVariantLists(geneData, geneList) {
    var geneVariants = {};
    var geneIDList = this.getGeneIDs(geneList);
    // initialize allGeneVariantsOrganized dict key value structure
    for (var x = 0; x < geneIDList.length; x++) {
      geneVariants[geneIDList[x]] = [];
    }
    for(var y = 0; y < geneData.length; y++) {
      geneVariants[geneData[y]["gene_id"]].push(geneData[y]["rsnum"]);
    }
    this.allGeneVariantsOrganized = geneVariants;
  }

  getGeneIDs(geneList) {
    var geneIDs = [];
    for (var i = 0; i < geneList.length; i++) {
      geneIDs.push(geneList[i]['gene_id']);
    }
    return geneIDs;
  }

  populatePopulationDropdown() {
    var populations = [
      {
        namecode: "AFR",
        name: "African",
        subPopulations: [
          { value: "YRI", viewValue: "Yoruba in Ibadan, Nigera" },
          { value: "LWK", viewValue: "Luhya in Webuye, Kenya" },
          { value: "GWD", viewValue: "Gambian in Western Gambia" },
          { value: "MSL", viewValue: "Mende in Sierra Leone" },
          { value: "ESN", viewValue: "Esan in Nigera" },
          { value: "ASW", viewValue: "Americans of African Ancestry in SW USA" },
          { value: "ACB", viewValue: "African Carribbeans in Barbados" },
        ]
      },
      {
        namecode: 'AMR',
        name: "Ad Mixed American",
        subPopulations: [
          { value: "MXL", viewValue: "Mexican Ancestry from Los Angeles, USA" },
          { value: "PUR", viewValue: "Puerto Ricans from Puerto Rico" },
          { value: "CLM", viewValue: "Colombians from Medellin, Colombia" },
          { value: "PEL", viewValue: "Peruvians from Lima, Peru" },
        ]
      },
      {
        namecode: "EAS",
        name: "East Asian",
        subPopulations: [
          { value: "CHB", viewValue: "Han Chinese in Bejing, China" },
          { value: "JPT", viewValue: "Japanese in Tokyo, Japan" },
          { value: "CHS", viewValue: "Southern Han Chinese" },
          { value: "CDX", viewValue: "Chinese Dai in Xishuangbanna, China" },
          { value: "KHV", viewValue: "Kinh in Ho Chi Minh City, Vietnam" },
        ]
      },
      {
        namecode: "EUR",
        name: "European",
        subPopulations: [
          { value: "CEU", viewValue: "Utah Residents from North and West Europe" },
          { value: "TSI", viewValue: "Toscani in Italia" },
          { value: "FIN", viewValue: "Finnish in Finland" },
          { value: "GBR", viewValue: "British in England and Scotland" },
          { value: "IBS", viewValue: "Iberian population in Spain" },
        ]
      },
      {
        namecode: "SAS",
        name: "South Asian",
        subPopulations: [
          { value: "GIH", viewValue: "Gujarati Indian from Houston, Texas" },
          { value: "PJL", viewValue: "Punjabi from Lahore, Pakistan" },
          { value: "BEB", viewValue: "Bengali from Bangladesh" },
          { value: "STU", viewValue: "Sri Lankan Tamil from the UK" },
          { value: "ITU", viewValue: "Indian Telugu from the UK" },
        ]
      }
    ];
    return populations;
  }

  selectAll() {
    if (this.selectedPop.length == 26) {
      this.selectedPop = [];
      if (this.populationSelectedAll == true) {
        this.populationSelectedAll = false;
      } else {
        this.populationSelectedAll = true;
      } 
    } else if (this.selectedPop.length < 26) {
      this.selectedPop = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"];
      if (this.populationSelectedAll == true) {
        this.populationSelectedAll = false;
      } else {
        this.populationSelectedAll = true;
      } 
    } else {
      // do nothing
    }
    this.inputChanged = true;
    this.recalculatePopAttempt = "true";
    this.returnPopulationGroupFinal();
  }

  unique(value, index, self) {
    return self.indexOf(value) === index;
  }

  containsAll(subarr, arr) {
    for (var i = 0, len = subarr.length; i < len; i++) {
      if (!arr.includes(subarr[i])) {
        return false;
      }
    }
    return true;
  }

  remove(element, src) {
    var newArray = JSON.parse(JSON.stringify(src));
    // console.log(newArray);
    for (var i = 0; i < newArray.length; i++) {
      var idx = -1;
      if (newArray[i] == element) {
        idx = i;
      }
      if (idx != -1) {
        newArray.splice(idx, 1);
      }
    }
    return newArray;
  }

  removeAll(subpop, src) {
    var newArray = JSON.parse(JSON.stringify(src));
    for (var i = 0; i < subpop.length; i++) {
      newArray = this.remove(subpop[i], newArray);
    }
    return newArray;
  }

  selectPopulationGroup(groupName) {
    var african = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"];
    var mixedAmerican = ["MXL", "PUR", "CLM", "PEL"];
    var eastAsian = ["CHB", "JPT", "CHS", "CDX", "KHV"];
    var european = ["CEU", "TSI", "FIN", "GBR", "IBS"];
    var southAsian = ["GIH", "PJL", "BEB", "STU", "ITU"];
    if (groupName == "AFR") {
      if (this.containsAll(african, this.selectedPop)) {
        this.selectedPop = this.removeAll(african, this.selectedPop);
        this.changePop();
      } else {
        this.selectedPop = (this.selectedPop.concat(african)).filter(this.unique);
        this.changePop();
      }
    }
    if (groupName == "AMR") {
      if (this.containsAll(mixedAmerican, this.selectedPop)) {
        this.selectedPop = this.removeAll(mixedAmerican, this.selectedPop);
        this.changePop();
      } else {
        this.selectedPop = (this.selectedPop.concat(mixedAmerican)).filter(this.unique);
        this.changePop();
      }
    }
    if (groupName == "EAS") {
      if (this.containsAll(eastAsian, this.selectedPop)) {
        this.selectedPop = this.removeAll(eastAsian, this.selectedPop);
        this.changePop();
      } else {
        this.selectedPop = (this.selectedPop.concat(eastAsian)).filter(this.unique);
        this.changePop();
      }
    }
    if (groupName == "EUR") {
      if (this.containsAll(european, this.selectedPop)) {
        this.selectedPop = this.removeAll(european, this.selectedPop);
        this.changePop();
      } else {
        this.selectedPop = (this.selectedPop.concat(european)).filter(this.unique);
        this.changePop();
      }
    }
    if (groupName == "SAS") {
      if (this.containsAll(southAsian, this.selectedPop)) {
        this.selectedPop = this.removeAll(southAsian, this.selectedPop);
        this.changePop();
      } else {
        this.selectedPop = (this.selectedPop.concat(southAsian)).filter(this.unique);
        this.changePop();
      }
    }
  }

  changePop() {
    if (this.selectedPop.length < 26) {
      this.populationSelectedAll = false;
    } else {
      this.populationSelectedAll = true;
    }
    this.inputChanged = true;
    this.recalculatePopAttempt = "true";
    this.returnPopulationGroupFinal();
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
        hoverData.push('chr' + geneData[i]['variant_id'] + '<br>' + geneData[i]['rsnum'] + '<br>' + 'Ref/Alt: ' + geneData[i]['ref'] + '/' + geneData[i]['alt'] + '<br>' + 'P-value: ' + geneData[i]['pval_nominal'] + '<br>' + 'Slope: ' + geneData[i]['slope']);
      } else {
        hoverData.push('chr' + geneData[i]['variant_id'] + '<br>' + 'Ref/Alt: ' + geneData[i]['ref'] + '/' + geneData[i]['alt'] + '<br>' + 'P-value: ' + geneData[i]['pval_nominal'] + '<br>' + 'Slope: ' + geneData[i]['slope']);
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

  getColorData(geneData) {
    var colorData = [];
    for (var i = 0; i < geneData.length; i++) {
      if (geneData[i]['R2']) {
        colorData.push(geneData[i]['R2']);
      } else {
        colorData.push(0.0);
      }
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

  locuszoomPlotGWAS(geneData, geneGWASData, geneDataRC, qDataTopAnnot) {
    var xData = this.getXData(geneData);
    var yData = this.getYData(geneData);
    var yGWASData = this.getYGWASData(geneGWASData);
    var colorData = this.getColorData(geneData);
    var xDataRC = this.getXDataRC(geneDataRC);
    var yDataRC = this.getYDataRC(geneDataRC);
    var hoverData = this.getHoverData(geneData);
    var hoverDataRC = this.getHoverDataRC(geneDataRC);

    // highlight top point
    var topAnnotHighlight = {
      x: [qDataTopAnnot['pos'] / 1000000.0],
      y: [Math.log10(qDataTopAnnot['pval_nominal']) * -1.0],
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 15,
        color: "red"
      },
      yaxis: 'y2'
    };

    // graph GWAS scatter
    var trace1 = {
      x: xData,
      y: yGWASData,
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
      // xaxis: 'x',
      yaxis: 'y'
    };

    // graph recombination rate line
    var trace2 = {
      x: xDataRC,
      y: yDataRC,
      text: hoverDataRC,
      hoverinfo: 'text',
      yaxis: 'y3',
      type: 'scatter',
      line: {
        color: 'blue',
        width: 1
      }
    };

    // graph scatter
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
      // xaxis: 'x',
      yaxis: 'y2'
    };

    // graph recombination rate line
    var trace4 = {
      x: xDataRC,
      y: yDataRC,
      text: hoverDataRC,
      hoverinfo: 'text',
      yaxis: 'y4',
      type: 'scatter',
      line: {
        color: 'blue',
        width: 1
      }
    };

    var pdata = [topAnnotHighlight, trace1, trace2, trace3, trace4];

    // round most significant pval to next whole number
    // var maxY = Math.ceil(Math.log10(qDataTopAnnot['pval_nominal']) * -1.0);
    var chromosome = qDataTopAnnot['chr'];

    var playout = {
      title: {
        text: 'Locuszoom Plot',
        xref: 'paper'
      },
      width: 1000,
      height: 1100,
      yaxis: {
        // range: [0, maxY],
        autorange: true,
        title: "GWAS -log10(P-value)",
        domain: [0, 0.48],
        zeroline: false
      },
      yaxis2: {
        // range: [0, maxY],
        autorange: true,
        title: "eQTL -log10(P-value)",
        domain: [0.52, 1],
        zeroline: false
      },
      yaxis3: {
        // range: [0, maxY * 10],
        autorange: true,
        title: 'GWAS Recombination Rate (cM/Mb)',
        titlefont: {
          color: 'blue'
        },
        tickfont: {
          color: 'blue'
        },
        overlaying: 'y',
        side: 'right',
        showgrid: false,
        dtick: 50,
        zeroline: false
      },
      yaxis4: {
        // range: [0, maxY * 10],
        autorange: true,
        title: 'eQTL Recombination Rate (cM/Mb)',
        titlefont: {
          color: 'blue'
        },
        tickfont: {
          color: 'blue'
        },
        overlaying: 'y2',
        side: 'right',
        showgrid: false,
        dtick: 50,
        zeroline: false
      },
      xaxis: {
        autorange: true,
        title: "Chromosome " + chromosome + " (Mb)",
        zeroline: false
      },
      images: [
        {
          x: 0,
          y: 1,
          sizex: 1.0,
          sizey: 1.0,
          source: environment.endpoint + "assets/images/eqtl_locuszoom_r2_legend.png",
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
        t: 140
      },
      showlegend: false,
      clickmode: 'event',
      hovermode: 'closest',
      // paper_bgcolor: "#D3D3D3"
    };
    
    return {
      data: pdata,
      layout: playout, 
      // divId: "eqtl-locuszoom-plot",
      config: {
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"]
      }
    };
  }

  locuszoomPlot(geneData, geneDataRC, qDataTopAnnot) {
    var xData = this.getXData(geneData);
    var yData = this.getYData(geneData);
    var colorData = this.getColorData(geneData);
    var xDataRC = this.getXDataRC(geneDataRC);
    var yDataRC = this.getYDataRC(geneDataRC);
    var hoverData = this.getHoverData(geneData);
    var hoverDataRC = this.getHoverDataRC(geneDataRC);

    // highlight top point
    var topAnnotHighlight = {
      x: [qDataTopAnnot['pos'] / 1000000.0],
      y: [Math.log10(qDataTopAnnot['pval_nominal']) * -1.0],
      hoverinfo: 'none',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 15,
        color: "red"
      }
    };

    // graph scatter
    var trace1 = {
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
      }
    };
    
    // graph recombination rate line
    var trace2 = {
      x: xDataRC,
      y: yDataRC,
      text: hoverDataRC,
      hoverinfo: 'text',
      yaxis: 'y2',
      type: 'scatter',
      line: {
        color: 'blue',
        width: 1
      }
    };
    
    var pdata = [topAnnotHighlight, trace1, trace2];

    // round most significant pval to next whole number
    // var maxY = Math.ceil(Math.log10(qDataTopAnnot['pval_nominal']) * -1.0);
    var chromosome = qDataTopAnnot['chr'];
    
    var playout = {
      title: {
        text: 'Locuszoom Plot',
        xref: 'paper'
      },
      width: 1000,
      height: 700,
      yaxis: {
        // range: [0, maxY],
        autorange: true,
        title: "eQTL -log10(P-value)",
        zeroline: false
      },
      yaxis2: {
        // range: [0, maxY * 10],
        autorange: true,
        title: 'eQTL Recombination Rate (cM/Mb)',
        titlefont: {
          color: 'blue'
        },
        tickfont: {
          color: 'blue'
        },
        overlaying: 'y',
        side: 'right',
        showgrid: false,
        dtick: 50,
        zeroline: false
      },
      xaxis: {
        autorange: true,
        title: "Chromosome " + chromosome + " (Mb)",
        zeroline: false
      },
      images: [
        {
          x: 0,
          y: 1,
          sizex: 1.0,
          sizey: 1.0,
          source: environment.endpoint + "assets/images/eqtl_locuszoom_r2_legend.png",
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
        t: 140
      },
      showlegend: false,
      clickmode: 'event',
      hovermode: 'closest',
      // annotations: []
      // annotations: [{
      //   x: 152.210492,
      //   y: 14.477790583737521,
      //   text: '<b>1:152210492</b>' + 
      //                 '<br>P Value: <b>3.3282e-15</b>' + 
      //                 '<br>Ref. Allele: <b> ATT</b>' + 
      //                 '<br>────────────────<br>' + 
      //                 '<a href="https://www.google.com"><b>Make LD Reference</b></a>' + 
      //         '<br><a href="https://www.google.com"><b>Show Boxplots</b></a>',
      //   align: "left",
      //   showarrow: true,
      //   clicktoshow: 'onout',
      //   visible: false,
      //   bordercolor: "black",
      //   bgcolor: "white",
      //   borderpad: 10, 
      //   ax: 90,
      //   ay: 0
      // }]
    };
    return {
      data: pdata,
      layout: playout,
      // divId: "eqtl-locuszoom-plot",
      config: {
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"]
      }
    };
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
    // console.log("Recalculate!");
    // console.log(boxplotData.point_index);
    // var selectedRefString = this.locuszoomData[boxplotData["point_index"]]["rsnum"];
    var selectedRefString = this.popoverData["rsnum"];
    // console.log(selectedRefString);
    var selectedGeneString = this.selectedGene;
    var selectedPopString = this.selectedPop.join('+');
    var recalculateAttempt = "true";
    var recalculatePop = this.recalculatePopAttempt;
    var recalculateGene = this.recalculateGeneAttempt;
    var recalculateRef = "true";
    this.inputChanged = false;
    // reset
    this.data.changeMainData(null);
    this.data.changeSelectedTab(0);
    // calculate
    this.data.recalculateMain(this.associationFile, this.expressionFile, this.genotypeFile, this.gwasFile, this.requestID, selectedPopString, selectedGeneString, selectedRefString, recalculateAttempt, recalculatePop, recalculateGene, recalculateRef)
      .subscribe(
        res => this.data.changeMainData(res),
        error => this.handleError(error)
      )
  }

  closeWarning() {
    this.warningMessage = "";
  }

  linkLDpop() {
    var selectedRefString = this.popoverData["rsnum"];
    var QTopAnnotRef = this.locuszoomDataQTopAnnot["rsnum"];
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

  linkGenomADBrowser() {
    var variant_id = this.popoverData["variant_id"].split(":");
    var chromosome = variant_id[0];
    var position = variant_id[1];
    var ref = this.popoverData["ref"];
    var alt = this.popoverData["alt"];
    var url = "http://gnomad.broadinstitute.org/variant/" + chromosome + "-" + position + "-" + ref + "-" + alt
    var win = window.open(url, '_blank');
    win.focus();
  }

  async showBoxplot() {
    if (this.popoverData) {
      // console.log("LOCUSZOOM BOXPLOT DIALOG MODULE OPENED");
      this.dialog.open(EqtlResultsLocuszoomBoxplotsComponent, {
        data: this.popoverData
      });
      this.closePopover();
      this.popoverData = null;
    }
  }

  // plotlyClick(event, plot: PlotComponent) {
  //   if (event.points) {
  //     var point = event.points[0];
  //     this.popoverData = this.locuszoomData[event.points[0].pointIndex];
  //     var newAnnotation = {
  //       x: point.xaxis.d2l(point.x),
  //       y: point.yaxis.d2l(point.y),
  //       text: '<b>' + this.popoverData['variant_id'] + '</b>' + 
  //             '<br>P Value: <b>' + this.popoverData['pval_nominal'] + '</b>' + 
  //             '<br>Ref. Allele: <b>' + this.popoverData['ref'] + '</b>' + 
  //             '<br>────────────────<br>' + 
  //             '<a (click)="makeLDRef()"><b>Make LD Reference</b></a>' + 
  //             '<br><a href="https://www.google.com"><b>Show Boxplots</b></a>',
  //       align: "left",
  //       showarrow: true,
  //       clicktoshow: 'onout',
  //       visible: false,
  //       bordercolor: "black",
  //       bgcolor: "white",
  //       borderpad: 10, 
  //       ax: 90,
  //       ay: 0,
  //     };
  //     console.log(newAnnotation);

  //     // plot.relayout('div', 'annotations[' + this.annotationIndex + ']', newAnnotation);
  //     console.log(this.graph);
  //     if (this.graph.layout.annotations) {
  //       this.graph.layout.annotations.push(newAnnotation);
  //     } else {
  //       this.graph.layout.annotations = [newAnnotation];
  //     }
  //   }
  // }

  clickPoint(event, plot: PlotComponent) {
    if (event.points) {
      // console.log(event);
      if (event.points[0].hasOwnProperty("marker.color")) {
        // console.log("SHOW MARKER");
        // only show popovers for scatter points not recomb line (points w/ markers)
        var top = event.event.pointerY;
        var left = event.event.pointerX;
        // console.log("event", event);
        // console.log("event.event", event.event);
        // console.log("pointerY", top);
        // console.log("pointerX", left);
        // console.log("pageY", event.event.pageY);
        // console.log("pageX", event.event.pageX);
        // console.log(event.points[0]);
        this.popoverData = this.locuszoomData[event.points[0].pointIndex];
        $('.popover').show();
        if (this.collapseInput) {
          // console.log("INPUT PANEL COLLAPSED");
          // GWAS scatter shown
          if (this.GWASData[0] && this.locuszoomScatterData[0]) {
            $('.popover').show().css({
              position: "absolute",
              top: top + 700, 
              left: left + 190
            });
          } else {
            $('.popover').show().css({
              position: "absolute",
              top: top, 
              left: left + 190
            });
          }
        } else {
          // console.log("INPUT PANEL SHOWN");
          // GWAS scatter shown
          if (this.GWASData[0] && this.locuszoomScatterData[0]) {
            $('.popover').show().css({
              position: "absolute",
              top: top + 700, 
              left: left + 25
            });
          } else {
            $('.popover').show().css({
              position: "absolute",
              top: top, 
              left: left + 25
            });
          }
        }
        this.showPopover = true;
      }
    } else {
      this.closePopover2(event);
    }
  }

  // hoverPoint(event, plot: PlotComponent) {
  //   console.log(event.points[0]);
  //   // Plotly.Fx.hover([{curveNumber:event.points[0].curveNumber - 2, pointNumber: event.points[0].pointIndex}]);
  // }

  

  handleError(error) {
    console.log(error);
    var errorTrimmed = error.error.trim().split('\n');
    // var errorMessage = errorTrimmed.slice(1, errorTrimmed.length - 1).join(' ');
    var errorMessage = errorTrimmed[2];
    console.log(errorMessage);
    this.data.changeErrorMessage(errorMessage);
  }

  returnPopulationGroupFinal() {
    this.selectedPopFinal = this.selectedPop;
    var african = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"];
    var mixedAmerican = ["MXL", "PUR", "CLM", "PEL"];
    var eastAsian = ["CHB", "JPT", "CHS", "CDX", "KHV"];
    var european = ["CEU", "TSI", "FIN", "GBR", "IBS"];
    var southAsian = ["GIH", "PJL", "BEB", "STU", "ITU"];
    if (this.containsAll(african, this.selectedPop)) {
      this.selectedPopFinal = this.removeAll(african, this.selectedPopFinal);
      this.selectedPopFinal.push("AFR");
    }
    if (this.containsAll(mixedAmerican, this.selectedPop)) {
      this.selectedPopFinal = this.removeAll(mixedAmerican, this.selectedPopFinal);
      this.selectedPopFinal.push("AMR");
    }
    if (this.containsAll(eastAsian, this.selectedPop)) {
      this.selectedPopFinal = this.removeAll(eastAsian, this.selectedPopFinal);
      this.selectedPopFinal.push("EAS");
    }
    if (this.containsAll(european, this.selectedPop)) {
      this.selectedPopFinal = this.removeAll(european, this.selectedPopFinal);
      this.selectedPopFinal.push("EUR");
    }
    if (this.containsAll(southAsian, this.selectedPop)) {
      this.selectedPopFinal = this.removeAll(southAsian, this.selectedPopFinal);
      this.selectedPopFinal.push("SAS");
    }
    return this.selectedPopFinal;
  }

  // async makeLDRefSearch() {
  //   if (this.geneVariantList.includes(this.rsnumSearch)) {
  //     var selectedRefString = this.rsnumSearch;
  //     // console.log("makeLDRefSearch()", selectedRefString);
  //     var selectedGeneString = this.selectedGene;
  //     var selectedPopString = this.selectedPop.join('+');
  //     var recalculateAttempt = "true";
  //     var recalculatePop = this.recalculatePopAttempt;
  //     var recalculateGene = this.recalculateGeneAttempt;
  //     var recalculateRef = "true";
  //     this.inputChanged = false;
  //     // reset
  //     this.data.changeMainData(null);
  //     this.data.changeSelectedTab(0);
  //     // calculate
  //     this.data.recalculateMain(this.associationFile, this.expressionFile, this.genotypeFile, this.gwasFile, this.requestID, selectedPopString, selectedGeneString, selectedRefString, recalculateAttempt, recalculatePop, recalculateGene, recalculateRef)
  //       .subscribe(
  //         res => this.data.changeMainData(res),
  //         error => this.handleError(error)
  //       )
  //   } else {
  //     this.warningMessage = this.rsnumSearch + " not found in Locuszoom gene data."
  //   }
  // }

  refGeneChange() {
    this.inputChanged = true;
    this.recalculateGeneAttempt = "true";
    for(var i = 0; i < this.topGeneVariants.length; i++) {
      if (this.topGeneVariants[i]['gene_id'] == this.selectedGene) {
        this.rsnumSearch = this.topGeneVariants[i]['rsnum'];
      }
    }
  }

  enableSearch(rsnumSearchInputValue) {
      this.inputChanged = true;
      this.recalculateRefAttempt = "true";
      this.rsnumSearch = rsnumSearchInputValue;
  }

  async recalculatePopGeneRef() {
    var selectedGeneString = this.selectedGene;
    var selectedPopString = this.selectedPop.join('+');
    var recalculateAttempt = "true";
    var recalculatePop = this.recalculatePopAttempt;
    var recalculateGene = this.recalculateGeneAttempt;
    this.inputChanged = false;
    if (this.allGeneVariantsOrganized[selectedGeneString].includes(this.rsnumSearch) || this.rsnumSearch.length == 0) {
      var selectedRefString = this.rsnumSearch;
      var recalculateRef = this.recalculateRefAttempt;
      if (this.rsnumSearch.length == 0) {
        selectedRefString = "false";
        recalculateRef = "false";
      }
      
      // reset
      this.data.changeMainData(null);
      this.data.changeSelectedTab(0);
      // calculate
      this.data.recalculateMain(this.associationFile, this.expressionFile, this.genotypeFile, this.gwasFile, this.requestID, selectedPopString, selectedGeneString, selectedRefString, recalculateAttempt, recalculatePop, recalculateGene, recalculateRef)
        .subscribe(
          res => this.data.changeMainData(res),
          error => this.handleError(error)
        );
      this.recalculateRefAttempt = "false";
      this.recalculatePopAttempt = "false";
      this.recalculateGeneAttempt = "false";
    } else {
      this.warningMessage = this.rsnumSearch + " not found in the association data file for the chosen reference gene. Please enter another variant."
    }
  }

  getScatterX(scatterData) {
    var p_values = [];
    for (var i = 0; i < scatterData.length; i++) {
      p_values.push(Math.log10(scatterData[i]['pvalue']) * -1.0);
    }
    return p_values;
  }

  getScatterY(scatterData) {
    var pval_nominals = [];
    for (var i = 0; i < scatterData.length; i++) {
      pval_nominals.push(Math.log10(scatterData[i]['pval_nominal']) * -1.0);
    }
    return pval_nominals;
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

  getLinearRegression(xData, yData) {
    // console.log("xData", xData);
    // console.log("yData", yData);
    // var sum = (accumulator, currentValue) => accumulator + currentValue;
    var n = (xData.length + yData.length) / 2;
    var xy = []
    var x2 = xData.map(x => x * x);
    for (var i = 0; i < xData.length; i++) {
      xy.push(xData[i] * yData[i]);
    }
    var sum_x = this.getSum(xData);
    // console.log("sum_x", sum_x);
    var sum_y = this.getSum(yData);
    // console.log("sum_y", sum_y);
    var sum_xy = this.getSum(xy);
    // console.log("sum_xy", sum_xy);
    var sum_x2 = this.getSum(x2);
    // console.log("sum_x2", sum_x2);
    var a_numer = (sum_y * sum_x2) - (sum_x * sum_xy);
    var a_denom = (n * sum_x2) - (sum_x * sum_x);
    var a = a_numer / a_denom;
    var b_numer = (n * sum_xy) - (sum_x * sum_y);
    var b_denom = (n * sum_x2) - (sum_x * sum_x);
    var b = b_numer / b_denom;
    return [a, b];
  }

  locuszoomScatterPlot(scatterData, scatterTitle) {
    var xData = this.getScatterX(scatterData);
    // console.log("xData", xData);
    var yData = this.getScatterY(scatterData);
    // console.log("yData", yData);
    var trace1 = {
      x: xData,
      y: yData,
      // hoverinfo: 'x+y',
      mode: 'markers',
      type: 'scatter',
      marker: {
        size: 5,
        color: "grey",
        line: {
          color: 'black',
          width: 1
        },
      }
    };
    // [a, b] -> Y = a + bX
    var linear_regression = this.getLinearRegression(xData, yData);
    var a = linear_regression[0];
    // console.log("a", a);
    var b = linear_regression[1];
    // console.log("b", b);
    var xMin = Math.min.apply(null, xData);
    // console.log("xMin", xMin);
    var xMax = Math.max.apply(null, xData);
    // console.log("xMax", xMax);
    var yMin = a + (b * xMin);
    // console.log("yMin", yMin);
    var yMax = a + (b * xMax);
    // console.log("yMax", yMax);
    var trace2 = {
      x: [xMin, xMax],
      y: [yMin, yMax],
      // hoverinfo: 'x+y',
      mode: 'lines',
      type: 'scatter',
      name: "lines",
      line: {
        color: "blue",
      }
    };
    var pdata = [trace1, trace2];
    // var pdata = [trace1];
    var playout = {
      title: {
        text: scatterTitle,
        xref: 'paper'
      },
      width: 1000,
      height: 700,
      yaxis: {
        autorange: true,
        title: "-log10(eQTL nominal pval)",
      },
      xaxis: {
        autorange: true,
        title: "-log10(GWAS pvalue)",
      },
      margin: {
        l: 40,
        r: 40,
        b: 80,
      },
      showlegend: false,
    };
    return {
      data: pdata,
      layout: playout,
      config: {
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "hoverCompareCartesian", "hoverClosestCartesian"]
      }
    };
  }

}

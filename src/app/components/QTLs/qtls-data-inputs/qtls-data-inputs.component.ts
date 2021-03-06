import { Component, OnInit, ChangeDetectorRef } from '@angular/core';
import { QTLsResultsService } from '../../../services/qtls-results.service';

import { FormGroup, FormControl, Validators } from '@angular/forms';

// declare var $; // declare jquery $
// import * as $ from 'jquery';
declare let $: any;

export interface GTExDataset {
  value: string;
  viewValue: string;
}

@Component({
  selector: 'app-qtls-data-inputs',
  templateUrl: './qtls-data-inputs.component.html',
  styleUrls: ['./qtls-data-inputs.component.css']
})

export class QTLsDataInputsComponent implements OnInit {

  qtlsForm = new FormGroup({
    associationFile: new FormControl({value: '', disabled: false}, Validators.required),
    quantificationFile: new FormControl({value: '', disabled: false}), 
    genotypeFile: new FormControl({value: '', disabled: false}), 
    gwasFile: new FormControl({value: '', disabled: false}),
    LDFile: new FormControl({value: '', disabled: false}),
    cisDistanceInput: new FormControl({value: "100", disabled: false}, [Validators.pattern("^(\-?(?!0)[0-9]+)?$"), Validators.min(1), Validators.max(2000), Validators.required]),
    rsnumber: new FormControl({value: '', disabled: false}, [Validators.pattern("^(rs[0-9]+)?$")])
  });

  mainData: Object;
  resultStatus: boolean;
  errorMessage: string;
  warningMessage: string;
  public resetColor = null;
  // selectLoadBoxplotData: boolean;
  // selectLoadGWASData: boolean;
  // selectLoadLDData: boolean;
  selectedTab: number;
  popoverData: Object;
  qtlsType: string;
  disableQTLsToggle: boolean;

  selectLoadQTLsSamples: boolean;
  selectLoadGWASSample: boolean;
  selectLoadLDSample: boolean;

  selectedDist: string;
  rsnumSearch: string;

  hyprcolocIsLoading: boolean;
  ecaviarIsLoading: boolean;

  GTExDatasets: GTExDataset[] = [
    {value: 'ds-0', viewValue: 'Dataset 1'},
    {value: 'ds-1', viewValue: 'Dataset 2'},
    {value: 'ds-2', viewValue: 'Dataset 3'}
  ];

  constructor(private cdr: ChangeDetectorRef, private data: QTLsResultsService) { }

  ngOnInit() {
    this.data.currentHyprcolocIsLoading.subscribe(hyprcolocIsLoading => {
      this.hyprcolocIsLoading = hyprcolocIsLoading
      if (!this.hyprcolocIsLoading) {
        $("#reset-tooltip").tooltip("enable");
        $("#reset-tooltip").tooltip("hide");
        $("#reset-tooltip").tooltip("disable");
      }
    });

    this.data.currentEcaviarIsLoading.subscribe(ecaviarIsLoading => {
      this.ecaviarIsLoading = ecaviarIsLoading
      if (!this.ecaviarIsLoading) {
        $("#reset-tooltip").tooltip("enable");
        $("#reset-tooltip").tooltip("hide");
        $("#reset-tooltip").tooltip("disable");
      }
    });

    this.data.currentMainData.subscribe(mainData => {
      this.mainData = mainData
      // populate LD Reference field with default variant after initial calculation if blank
      if (mainData && (this.rsnumSearch == null || this.rsnumSearch.length == 0)) {
        var locusAlignmentDataQTopAnnot = mainData["locus_alignment"]["top"][0][0]; // locus alignment Top Gene data
        this.rsnumSearch = locusAlignmentDataQTopAnnot["rsnum"];
        this.qtlsForm.value.rsnumber = locusAlignmentDataQTopAnnot["rsnum"];
      }
    });
    this.data.currentResultStatus.subscribe(resultStatus => this.resultStatus = resultStatus);
    this.data.currentSelectedTab.subscribe(selectedTab => this.selectedTab = selectedTab);
    this.data.currentQtlsType.subscribe(qtlsType => {
      this.qtlsType = qtlsType;
    });
    this.data.currentErrorMessage.subscribe(errorMessage => {
      this.errorMessage = errorMessage;
      if (this.errorMessage) {
        this.resetColor = 'warn';
      } else {
        this.resetColor = null;
      }
    });
    this.disableQTLsToggle = false;
    this.selectLoadQTLsSamples = false;
    this.selectLoadGWASSample = false;
    this.selectLoadLDSample = false;
    this.selectedDist = "100"; // default cis-QTL distance (in Kb)
  }

  toggleAssocGTEx() {
    console.log(this.qtlsType);
    this.clearAssociationFile();
  }

  // downloadQTLsSamples() {
  //   console.log("Download samples...");
  //   var sampleFiles = [
  //     "assets/files/MX2.eQTL.txt",
  //     "assets/files/MX2.quantification.txt",
  //     "assets/files/MX2.genotyping.txt",
  //     "assets/files/MX2.LD.gz"
  //   ];
  //   for (var i = 0; i < sampleFiles.length; i ++) {
  //     var a = document.createElement("a");
  //     a.href = sampleFiles[i];
  //     a.download = sampleFiles[i].split('\/').pop();
  //     a.click();
  //   }
  // }

  clearAssociationFile() {
    this.qtlsForm.setControl('associationFile', new FormControl({value: '', disabled: false}, Validators.required));
    this.qtlsForm.value.associationFile = false;
    $("#association-file").val("");
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
  }

  clearQuantificationFile() {
    this.qtlsForm.setControl('quantificationFile', new FormControl({value: '', disabled: false}));
    this.qtlsForm.value.quantificationFile = false;
    $("#quantification-file").val("");
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
  }

  clearGenotypeFile() {
    this.qtlsForm.setControl('genotypeFile', new FormControl({value: '', disabled: false}));
    this.qtlsForm.value.genotypeFile = false;
    $("#genotype-file").val("");
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
  }

  clearGWASFile() {
    this.qtlsForm.setControl('gwasFile', new FormControl({value: '', disabled: false}));
    this.qtlsForm.value.gwasFile = false;
    $("#gwas-file").val("");
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
  }

  clearLDFile() {
    this.qtlsForm.setControl('LDFile', new FormControl({value: '', disabled: false}));
    this.qtlsForm.value.LDFile = false;
    $("#LD-file").val("");
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
  }

  loadQTLsSampleDataFiles() { // if user unloads sample QTLs data files
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
    if (this.selectLoadQTLsSamples == true) {
      this.selectLoadQTLsSamples = false;
      this.disableQTLsToggle = false;
      this.qtlsForm.setControl('associationFile', new FormControl({value: '', disabled: false}, Validators.required));
      this.qtlsForm.setControl('quantificationFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('genotypeFile', new FormControl({value: '', disabled: false}));
      $("#qtls-data-input-association-file").removeClass("disabled-overlay");
      $("#qtls-data-input-quantification-file").removeClass("disabled-overlay");
      $("#qtls-data-input-genotype-file").removeClass("disabled-overlay");
      this.qtlsForm.value.associationFile = false;
      this.qtlsForm.value.quantificationFile = false;
      this.qtlsForm.value.genotypeFile = false;
      $("#association-file").val("");
      $("#quantification-file").val("");
      $("#genotype-file").val("");
      this.data.changeDisableLocusQuantification(true);
      this.loadLDSampleDataFile();
      // if (this.selectLoadGWASSample == true) {
      //   this.loadGWASSampleDataFile(); // toggle load GWAS data file
      // }
    } else { // if user loads sample QTLs data files
      this.selectLoadQTLsSamples = true;
      this.qtlsType = "assoc";
      this.disableQTLsToggle = true;
      this.qtlsForm.setControl('associationFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('quantificationFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('genotypeFile', new FormControl({value: '', disabled: false}));
      $("#qtls-data-input-association-file").addClass("disabled-overlay");
      $("#qtls-data-input-quantification-file").addClass("disabled-overlay");
      $("#qtls-data-input-genotype-file").addClass("disabled-overlay");
      this.qtlsForm.value.associationFile = false;
      this.qtlsForm.value.quantificationFile = false;
      this.qtlsForm.value.genotypeFile = false;
      $("#association-file").val("");
      $("#quantification-file").val("");
      $("#genotype-file").val("");
      this.data.changeDisableLocusQuantification(false);
      this.loadLDSampleDataFile();
      // if (this.selectLoadGWASSample == false) {
      //   this.loadGWASSampleDataFile(); // toggle load GWAS data file
      // }
    }
  }

  loadGWASSampleDataFile() { // if user unloads sample GWAS data file
    $("#quantification-file-tooltip").tooltip("enable");
    $("#quantification-file-tooltip").tooltip("hide");
    $("#quantification-file-tooltip").tooltip("disable");
    $("#genotype-file-tooltip").tooltip("enable");
    $("#genotype-file-tooltip").tooltip("hide");
    $("#genotype-file-tooltip").tooltip("disable");
    if (this.selectLoadGWASSample == true) {
      this.selectLoadGWASSample = false;
      this.qtlsForm.setControl('gwasFile', new FormControl({value: '', disabled: false}));
      $("#qtls-data-input-gwas-file").removeClass("disabled-overlay");
      this.qtlsForm.value.gwasFile = false;
      $("#gwas-file").val("");
      this.data.changeDisableLocusColocalization(true);
    } else { // if user loads sample GWAS data file
      this.selectLoadGWASSample = true;
      this.qtlsForm.setControl('gwasFile', new FormControl({value: '', disabled: false}));
      $("#qtls-data-input-gwas-file").addClass("disabled-overlay");
      this.qtlsForm.value.gwasFile = false;
      $("#gwas-file").val("");
      this.data.changeDisableLocusColocalization(false);
    }
  }

  loadLDSampleDataFile() {
    if (this.selectLoadLDSample == true) { // if unloads sample LD data file
      this.selectLoadLDSample = false;
      this.qtlsForm.setControl('LDFile', new FormControl({value: '', disabled: false}));
      $("#qtls-data-input-LD-file").removeClass("disabled-overlay");
      this.qtlsForm.value.LDFile = false;
      $("#LD-file").val("");
    } else { // if user loads sample LD data file
      this.selectLoadLDSample = true;
      this.qtlsForm.setControl('LDFile', new FormControl({value: '', disabled: false}));
      $("#qtls-data-input-LD-file").addClass("disabled-overlay");
      this.qtlsForm.value.LDFile = false;
      $("#LD-file").val("");
    }
  }

  handleError(error) {
    var errorMessage = "";
    if (error.error.includes("Rscript")) {
      var errorTrimmed = error.error.trim().split('\n');
      errorMessage = errorTrimmed[2];
    } else {
      errorMessage = error.error;
    }
    this.data.changeErrorMessage(errorMessage);
  }

  enableSearchCISDistance(event: any) {
    // this.inputChanged = true;
    // this.recalculateDistAttempt = "true";
    this.selectedDist = event.target.value;
    this.qtlsForm.value.cisDistanceInput = event.target.value;
  }

  enableSearchLDRef(event: any) {
    // this.inputChanged = true;
    // this.recalculateRefAttempt = "true";
    this.rsnumSearch = event.target.value;
    this.qtlsForm.value.rsnumber = event.target.value;
  }

  clearCISDistField() {
    this.selectedDist = null;
    this.qtlsForm.value.cisDistanceInput = '';
  }

  clearLDRefField() {
    this.rsnumSearch = '';
    this.qtlsForm.value.rsnumber = '';
  }

  cisDistErrorMsg() {
    var msg = "";
    if (this.qtlsForm.value.cisDistanceInput > 2000) {
      msg = "Distance must be <= 2000 Kb";
    } else if (this.qtlsForm.value.cisDistanceInput < 1) {
      msg = "Distance must be >= 1 Kb";
    } else {
      msg = "Invalid cis-QTL Distance";
    }
    if (this.qtlsForm.value.cisDistanceInput == null || this.qtlsForm.value.cisDistanceInput == '') {
      msg = "Input required";
    }
    return msg;
  }

  async submit() {
    var request_id = Date.now().toString();
    var selectedDistNumber = this.selectedDist;
    var { associationFile, quantificationFile, genotypeFile, gwasFile, LDFile } = this.qtlsForm.value;
    var formData = new FormData();
    // custom tooltip validation - if quantification file is submitted, need genotype file and vice versa. All or none.
    if ((quantificationFile == null || quantificationFile == false || quantificationFile.length == 0) && (genotypeFile != null && genotypeFile.length > 0)) {
      // show tooltip on quantification file if user tries to calculate with genotype file but no quantification file
      // prevent calculation
      // console.log("no quantification - yes genotype"); 
      $("#quantification-file-tooltip").tooltip("enable");
      $("#quantification-file-tooltip").tooltip("show");
      $("#quantification-file-tooltip").tooltip("disable");
    } else if ((quantificationFile != null && quantificationFile.length > 0) && (genotypeFile == null || genotypeFile == false || genotypeFile.length == 0)) {
      // show tooltip on genotype file if user tries to calculate with quantification file but no genotype file
      // prevent calculation
      // console.log("yes quantification - no genotype");
      $("#genotype-file-tooltip").tooltip("enable");
      $("#genotype-file-tooltip").tooltip("show");
      $("#genotype-file-tooltip").tooltip("disable");
    } else {
      // only run calculation if both or none quantification and genotype type files are inputted
      // console.log("yes/no quantification - yes/no genotype");
      $("#quantification-file-tooltip").tooltip("enable");
      $("#quantification-file-tooltip").tooltip("hide");
      $("#quantification-file-tooltip").tooltip("disable");
      $("#genotype-file-tooltip").tooltip("enable");
      $("#genotype-file-tooltip").tooltip("hide");
      $("#genotype-file-tooltip").tooltip("disable");
      this.data.changeResultStatus(true);
      this.data.changeBlurLoadMain(true);
      $(".blur-loading-main").addClass("blur-overlay");
      $(".disabled-post-calc").addClass("disabled-overlay");
      $("#qtls-data-input-association-file").addClass("disabled-overlay");
      $("#qtls-data-input-quantification-file").addClass("disabled-overlay");
      $("#qtls-data-input-genotype-file").addClass("disabled-overlay");
      $("#qtls-data-input-gwas-file").addClass("disabled-overlay");
      $("#qtls-data-input-LD-file").addClass("disabled-overlay");
      this.qtlsForm.controls['cisDistanceInput'].disable();
      this.qtlsForm.controls['rsnumber'].disable();
      formData.append('request_id', request_id); // generate calculation request ID
      // formData.append('request_id', Date.now().toString()); // generate calculation request ID
      formData.append('select_pop', "false"); // set default population to 'false' -> 'EUR' in R
      formData.append('select_gene', "false"); // set default gene to 'false' -> QData top gene in R

      formData.append('select_dist', selectedDistNumber); // required input field
      if (this.rsnumSearch == null || this.rsnumSearch.length == 0) {
        // optional input field
        // if empty, set default rsnum to 'false' -> QData top gene's rsnum in R
        formData.append('select_ref', "false");
      } else {
        formData.append('select_ref', this.rsnumSearch);
      }
      
      formData.append('recalculateAttempt', "false");
      formData.append('recalculatePop', "false");
      formData.append('recalculateGene', "false");
      formData.append('recalculateDist', "false");
      formData.append('recalculateRef', "false");
      if (this.selectLoadQTLsSamples) {
        formData.append('select_qtls_samples', "true");
        // var select_qtls_samples = "true";
      } else {
        formData.append('select_qtls_samples', "false");
        // var select_qtls_samples = "false";
      }
      if (this.selectLoadGWASSample) {
        formData.append('select_gwas_sample', "true");
        // var select_gwas_sample = "true";
      } else {
        formData.append('select_gwas_sample', "false");
        // var select_gwas_sample = "false";
      }

      if (associationFile != null && associationFile.length > 0) {
        formData.append('association-file', associationFile[0]);
      }

      if ((quantificationFile != null && quantificationFile.length > 0) && (genotypeFile != null && genotypeFile.length > 0)) {
        formData.append('quantification-file', quantificationFile[0]);
        formData.append('genotype-file', genotypeFile[0]);
        this.data.changeDisableLocusQuantification(false);
      } else if (this.selectLoadQTLsSamples) {
        this.data.changeDisableLocusQuantification(false);
      } else {
        this.data.changeDisableLocusQuantification(true);
      }
      
      if (gwasFile != null && gwasFile.length > 0) {
        formData.append('gwas-file', gwasFile[0]);
        this.data.changeDisableLocusColocalization(false);
      } else if (this.selectLoadGWASSample) {
        this.data.changeDisableLocusColocalization(false);
      } else {
        this.data.changeDisableLocusColocalization(true);
      }

      if (LDFile != null && LDFile.length > 0) {
        formData.append('LD-file', LDFile[0]);
      }

      // console.log("RUN MAIN CALCULATION");
      this.data.calculateMain(formData)
        .subscribe(
          res => {
            // console.log("RESPONSE MAIN");
            this.data.changeMainData(res);
            this.data.changeBlurLoadMain(false);
            $(".blur-loading-main").removeClass("blur-overlay");
            $("#qtls-data-input-association-file").removeClass("disabled-overlay");
            $("#qtls-data-input-quantification-file").removeClass("disabled-overlay");
            $("#qtls-data-input-genotype-file").removeClass("disabled-overlay");
            $("#qtls-data-input-gwas-file").removeClass("disabled-overlay");
            $("#qtls-data-input-LD-file").removeClass("disabled-overlay");
          },
          error => {
            this.handleError(error);
            this.data.changeBlurLoadMain(false);
            $(".blur-loading-main").removeClass("blur-overlay");
            $("#qtls-data-input-association-file").removeClass("disabled-overlay");
            $("#qtls-data-input-quantification-file").removeClass("disabled-overlay");
            $("#qtls-data-input-genotype-file").removeClass("disabled-overlay");
            $("#qtls-data-input-gwas-file").removeClass("disabled-overlay");
            $("#qtls-data-input-LD-file").removeClass("disabled-overlay");
          }
        );
    }
  } 

  reset() {
    if (this.hyprcolocIsLoading || this.ecaviarIsLoading) {
      // show reset tooltip warning and prevent reset if any hyprcoloc/ecaviar calculation is ongoing
      $("#reset-tooltip").tooltip("enable");
      $("#reset-tooltip").tooltip("show");
      $("#reset-tooltip").tooltip("disable");
    } else {
      // hide recalculate/reset tooltip & reset hyprcoloc/ecaviar loading status
      $("#reset-tooltip").tooltip("enable");
      $("#reset-tooltip").tooltip("hide");
      $("#reset-tooltip").tooltip("disable");
      $("#recalculate-tooltip").tooltip("enable");
      $("#recalculate-tooltip").tooltip("hide");
      $("#recalculate-tooltip").tooltip("disable");
      this.data.changeHyprcolocIsLoading(false);
      this.data.changeEcaviarIsLoading(false);
      // hide quantification file genotype file error tooltips if displayed
      $("#quantification-file-tooltip").tooltip("enable");
      $("#quantification-file-tooltip").tooltip("hide");
      $("#quantification-file-tooltip").tooltip("disable");
      $("#genotype-file-tooltip").tooltip("enable");
      $("#genotype-file-tooltip").tooltip("hide");
      $("#genotype-file-tooltip").tooltip("disable");
      this.disableQTLsToggle = false;
      this.data.changeResultStatus(false);
      this.data.changeBlurLoadMain(false);
      // this.data.changeBlurLoadECAVIAR(false);
      $(".blur-loading-main").removeClass("blur-overlay");
      $(".disabled-post-calc").removeClass("disabled-overlay");
      // remove all calculated data
      this.data.changeMainData(null);
      // remove all eCAVIAR calculated data
      this.data.changeECAVIARData(null);
      // remove all Hyprcoloc calculate data
      this.data.changeHyprcolocData(null);
      // remove any error messages
      this.data.changeErrorMessage('');
      // choose default association data file toggle
      // this.qtlsType = "assoc";
      this.data.changeQtlsType("assoc");
      // change result tab back to Locus Alignment
      this.data.changeSelectedTab(0);
      // clear any file inputs
      this.qtlsForm.setControl('associationFile', new FormControl({value: '', disabled: false}, Validators.required));
      this.qtlsForm.setControl('quantificationFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('genotypeFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('gwasFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('LDFile', new FormControl({value: '', disabled: false}));
      this.qtlsForm.setControl('cisDistanceInput', new FormControl({value: '100', disabled: false}, [Validators.pattern("^(\-?(?!0)[0-9]+)?$"), Validators.min(1), Validators.max(2000), Validators.required]));
      this.qtlsForm.setControl('rsnumber', new FormControl({value: '', disabled: false}, [Validators.pattern("^(rs[0-9]+)?$")]));
      this.qtlsForm.value.associationFile = false;
      this.qtlsForm.value.quantificationFile = false;
      this.qtlsForm.value.genotypeFile = false;
      this.qtlsForm.value.gwasFile = false;
      this.qtlsForm.value.LDFile = false;
      $("#association-file").val("");
      $("#quantification-file").val("");
      $("#genotype-file").val("");
      $("#gwas-file").val("");
      $("#LD-file").val("");
      // disable locus colocalization and locus quantification result tabs
      this.data.changeDisableLocusColocalization(true);
      this.data.changeDisableLocusQuantification(true);
      $("#qtls-data-input-association-file").removeClass("disabled-overlay");
      $("#qtls-data-input-quantification-file").removeClass("disabled-overlay");
      $("#qtls-data-input-genotype-file").removeClass("disabled-overlay");
      $("#qtls-data-input-gwas-file").removeClass("disabled-overlay");
      $("#qtls-data-input-LD-file").removeClass("disabled-overlay");
      // clear any loaded sample files
      this.selectLoadQTLsSamples = false;
      this.selectLoadGWASSample = false;
      this.selectLoadLDSample = false;
      this.selectedDist = "100";
      this.rsnumSearch = "";
    }
  }

}

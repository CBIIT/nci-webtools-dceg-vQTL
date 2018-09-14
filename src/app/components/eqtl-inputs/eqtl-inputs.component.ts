import { Component, OnInit } from '@angular/core';
import { EqtlResultsComponent } from '../eqtl-results/eqtl-results.component';
import { EqtlResultsService } from '../../services/eqtl-results.service';

import { FormGroup, FormControl } from '@angular/forms';
import { Validators } from '@angular/forms';

import { environment } from '../../../environments/environment';

@Component({
  selector: 'app-eqtl-inputs',
  templateUrl: './eqtl-inputs.component.html',
  styleUrls: ['./eqtl-inputs.component.css']
})

export class EqtlInputsComponent implements OnInit {

  eqtlForm = new FormGroup({
    expressionFile: new FormControl('', Validators.required), 
    genotypeFile: new FormControl('', Validators.required), 
    associationFile: new FormControl('', Validators.required)
  });

  message: Object;
  resultStatus: boolean;
  errorMessage: string;
  public resetColor = null;

  constructor(private data: EqtlResultsService) { }

  ngOnInit() {
    this.eqtlForm.valueChanges.subscribe(formValue => {
      console.log(formValue);
    });

    this.data.currentMessage.subscribe(message => this.message = message);
    this.data.currentResultStatus.subscribe(resultStatus => this.resultStatus = resultStatus);
    this.data.currentErrorMessage.subscribe(errorMessage => {
      this.errorMessage = errorMessage;
      if (this.errorMessage) {
        this.resetColor = 'warn';
      } else {
        this.resetColor = null;
      }
    });
  }

  // this.data.currentMessage.subscribe(message => {
  //   this.message = message;
  //   this.graph = this.exampleBoxPlot(this.message);
  // });

  // newMessage() {
  //   this.data.changeMessage("Hello from Sibling");
  // }

  async submit() {
    this.data.changeResultStatus(true);

    const { expressionFile, genotypeFile, associationFile } = this.eqtlForm.value;
    console.log([expressionFile[0].name, genotypeFile[0].name, associationFile[0].name]);

    const formData = new FormData();
    formData.append('expression-file', expressionFile[0]);
    formData.append('genotype-file', genotypeFile[0]);
    formData.append('association-file', associationFile[0]);

    this.data.getResults(formData)
      .subscribe(
        res => this.data.changeMessage(res),
        error => this.handleError(error)
      )
    
    // const response = await fetch(environment.endpoint + '/upload-file', {
    //   method: 'POST',
    //   body: formData
    // });

    // // console.log(await response.text());

    // this.data.changeMessage(await response.json());
  } 

  handleError(error) {
    var errorTrimmed = error.error.trim().split('\n');
    // var errorMessage = errorTrimmed.slice(1, errorTrimmed.length - 1).join(' ');
    var errorMessage = errorTrimmed[2];
    console.log(errorMessage);
    this.data.changeErrorMessage(errorMessage);
  }

  reset() {
    this.data.changeResultStatus(false);
    this.data.changeMessage('');
    this.data.changeErrorMessage('');
  }

  // onSubmit() {
  //   // TODO: Use EventEmitter with form value
  //   console.log("Form has been submitted.");
  //   console.warn(this.eqtlForm.value);
  // }
}
